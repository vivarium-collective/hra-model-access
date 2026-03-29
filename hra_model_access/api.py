"""Cached HTTP clients for BioModels, OLS, HGNC, UniProt, PubMed + SBML parsing."""

from __future__ import annotations

import json
import re
import shelve
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree as ET

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

BIOMODELS = "https://www.ebi.ac.uk/biomodels"
OLS = "https://www.ebi.ac.uk/ols4/api"
HGNC = "https://rest.genenames.org"
UNIPROT = "https://rest.uniprot.org/uniprotkb"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

_CACHE_DIR = Path.home() / ".cache" / "hra_model_access"
_session: requests.Session | None = None

SUBCELLULAR = frozenset({
    "cytosol", "cytoplasm", "nucleus", "mitochondria", "mitochondrion",
    "endoplasmic reticulum", "golgi", "golgi apparatus", "lysosome",
    "peroxisome", "plasma membrane", "membrane", "extracellular",
    "extracellular space", "default", "compartment", "cell",
})

_RDF_LI = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}li"
_RDF_RES = "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource"

_SBML_URI_RE = [
    ("uniprot", re.compile(r"identifiers\.org/uniprot/([A-Z0-9]+)")),
    ("go",      re.compile(r"identifiers\.org/go/(GO:\d+)")),
    ("chebi",   re.compile(r"identifiers\.org/chebi/(CHEBI:\d+)")),
    ("bto",     re.compile(r"identifiers\.org/bto/(BTO:\d+)")),
    ("cl",      re.compile(r"identifiers\.org/cl/(CL:\d+)")),
    ("doid",    re.compile(r"identifiers\.org/(?:doid|DOID)/(DOID:\d+)")),
]


# -- Session / cache ---------------------------------------------------------

def _session_get() -> requests.Session:
    global _session
    if _session is None:
        _session = requests.Session()
        _session.mount("https://", HTTPAdapter(
            max_retries=Retry(total=3, backoff_factor=0.5,
                              status_forcelist=[429, 500, 502, 503, 504])))
    return _session


def _cache() -> shelve.Shelf:
    if not hasattr(_cache, "_shelf"):
        _CACHE_DIR.mkdir(parents=True, exist_ok=True)
        _cache._shelf = shelve.open(str(_CACHE_DIR / "api_cache"))
    return _cache._shelf


def _get(url: str, *, params: dict | None = None, headers: dict | None = None,
         mode: str = "json"):
    """Cached GET. mode='json' returns dict|None, 'text' returns str|None, 'bytes' returns bytes|None."""
    key = f"{mode}|{url}|{json.dumps(params or {}, sort_keys=True)}"
    c = _cache()
    if key in c:
        return c[key]
    try:
        r = _session_get().get(url, params=params or {}, headers=headers or {},
                               timeout=20 if mode != "bytes" else 30)
        r.raise_for_status()
        val = r.json() if mode == "json" else (r.text if mode == "text" else r.content)
        c[key] = val
        return val
    except Exception:
        return None


# -- BioModels ---------------------------------------------------------------

def search_biomodels(query: str, max_results: int = 50) -> list[str]:
    ids, offset = [], 0
    while len(ids) < max_results:
        data = _get(f"{BIOMODELS}/search",
                    params={"query": query, "numResults": 10, "offset": offset, "format": "json"})
        if not data:
            break
        models = data.get("models", [])
        if not models:
            break
        ids.extend(m["id"] for m in models)
        offset += len(models)
        if offset >= data.get("matches", 0):
            break
    return ids[:max_results]


def fetch_model(model_id: str) -> dict:
    data = _get(f"{BIOMODELS}/{model_id}", params={"format": "json"})
    if data is None:
        raise ValueError(f"Failed to fetch {model_id}")
    return data


# -- SBML parsing ------------------------------------------------------------

@dataclass
class SBMLData:
    uniprot: set[str] = field(default_factory=set)
    go: set[str] = field(default_factory=set)
    chebi: set[str] = field(default_factory=set)
    bto: set[str] = field(default_factory=set)
    cl: set[str] = field(default_factory=set)
    doid: set[str] = field(default_factory=set)
    compartments: list[str] = field(default_factory=list)
    subcellular: list[str] = field(default_factory=list)


def fetch_sbml(details: dict) -> SBMLData:
    model_id = details.get("publicationId") or details.get("submissionId", "")
    sbml_files = [f for f in details.get("files", {}).get("main", [])
                  if f["name"].endswith((".xml", ".sbml"))]
    if not sbml_files or not model_id:
        return SBMLData()
    content = _get(f"{BIOMODELS}/model/download/{model_id}",
                   params={"filename": sbml_files[0]["name"]}, mode="bytes")
    if not content:
        return SBMLData()
    return _parse_sbml(content)


def _parse_sbml(raw: bytes) -> SBMLData:
    data = SBMLData()
    try:
        root = ET.fromstring(raw)
    except ET.ParseError:
        return data
    sets = {k: getattr(data, k) for k in ("uniprot", "go", "chebi", "bto", "cl", "doid")}
    for li in root.iter(_RDF_LI):
        uri = li.get(_RDF_RES, "")
        for key, pat in _SBML_URI_RE:
            m = pat.search(uri)
            if m:
                sets[key].add(m.group(1))
    for elem in root.iter():
        tag = elem.tag.rsplit("}", 1)[-1]
        if tag == "compartment":
            name = (elem.get("name") or elem.get("id") or "").strip()
            if not name or len(name) <= 2:
                continue
            (data.subcellular if name.lower() in SUBCELLULAR else data.compartments).append(name)
    return data


# -- PubMed ------------------------------------------------------------------

def fetch_abstract(pubmed_id: str) -> str:
    if not pubmed_id:
        return ""
    return _get(f"{EUTILS}/efetch.fcgi",
                params={"db": "pubmed", "id": pubmed_id, "rettype": "abstract",
                        "retmode": "text"},
                mode="text") or ""


# -- OLS ---------------------------------------------------------------------

def ols_lookup(term: str, ontology: str, prefix: str) -> tuple[str, str]:
    """Search OLS → (label, obo_id) or (term, "")."""
    data = _get(f"{OLS}/search", params={"q": term, "ontology": ontology, "rows": 3})
    for doc in (data or {}).get("response", {}).get("docs", []):
        obo_id = doc.get("obo_id", "")
        if obo_id.startswith(prefix):
            return doc.get("label", term), obo_id
    return term, ""


def lookup_uberon(name: str) -> tuple[str, str]:
    return ols_lookup(name, "uberon", "UBERON:")


def lookup_cl(name: str) -> tuple[str, str]:
    return ols_lookup(name, "cl", "CL:")


# -- HGNC / UniProt ----------------------------------------------------------

def lookup_hgnc(symbol: str) -> tuple[str, str]:
    data = _get(f"{HGNC}/fetch/symbol/{symbol}", headers={"Accept": "application/json"})
    for doc in (data or {}).get("response", {}).get("docs", []):
        return doc.get("symbol", symbol), doc.get("hgnc_id", "")
    return symbol, ""


def lookup_uniprot(accession: str) -> dict:
    data = _get(f"{UNIPROT}/{accession}.json", headers={"Accept": "application/json"})
    if not data:
        return {}
    result = {}
    genes = data.get("genes", [])
    if genes:
        sym = genes[0].get("geneName", {}).get("value", "")
        if sym:
            result["gene_symbol"] = sym
    for xref in data.get("uniProtKBCrossReferences", []):
        if xref.get("database") == "HGNC":
            result["hgnc_id"] = xref.get("id", "")
            break
    for c in data.get("comments", []):
        if c.get("commentType") == "TISSUE SPECIFICITY":
            texts = c.get("texts", [])
            if texts:
                result["tissue_specificity"] = texts[0].get("value", "")
    return result
