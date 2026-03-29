"""Enrich BioModels metadata with HRA ontology mappings.

Data sources: SBML annotations, SBML compartments, model-level annotations,
PubMed abstracts, UniProt tissue specificity, description keywords, disease inference.
"""

from __future__ import annotations

import re

from . import api

# -- Keyword tables ----------------------------------------------------------

ANATOMY_TERMS: dict[str, str] = {
    "liver": "liver", "hepatic": "liver",
    "pancreas": "pancreas", "pancreatic": "pancreas",
    "islet": "islet of Langerhans",
    "kidney": "kidney", "renal": "kidney",
    "heart": "heart", "cardiac": "heart",
    "lung": "lung", "pulmonary": "lung",
    "brain": "brain",
    "adipose": "adipose tissue",
    "skeletal muscle": "skeletal muscle tissue", "muscle": "skeletal muscle tissue",
    "blood": "blood",
    "bone marrow": "bone marrow", "spleen": "spleen", "thymus": "thymus",
    "intestin": "intestine", "colon": "colon",
    "stomach": "stomach", "gastric": "stomach",
    "skin": "skin",
    "breast": "breast", "mammary": "mammary gland",
    "prostate": "prostate gland", "ovary": "ovary", "uterus": "uterus",
    "testis": "testis", "thyroid": "thyroid gland",
    "adrenal": "adrenal gland", "pituitary": "pituitary gland",
    "retina": "retina",
    "tumor": "neoplasm", "tumour": "neoplasm",
}

CELL_TERMS: dict[str, str] = {
    "hepatocyte": "hepatocyte", "adipocyte": "adipocyte",
    "myocyte": "myocyte", "cardiomyocyte": "cardiomyocyte",
    "neuron": "neuron", "astrocyte": "astrocyte",
    "beta cell": "pancreatic beta cell", "beta-cell": "pancreatic beta cell",
    "macrophage": "macrophage", "neutrophil": "neutrophil",
    "fibroblast": "fibroblast",
    "endothelial cell": "endothelial cell", "epithelial cell": "epithelial cell",
    "erythrocyte": "erythrocyte",
    "osteoblast": "osteoblast", "osteoclast": "osteoclast",
    "keratinocyte": "keratinocyte", "melanocyte": "melanocyte",
}

# keyword → (context_label, associated_organs)
# Order matters: more specific patterns first
DISEASE_CONTEXT: list[tuple[str, str, list[str]]] = [
    ("type 2 diabetes",    "type 2 diabetes",    ["pancreas", "adipose tissue", "liver"]),
    ("type 1 diabetes",    "type 1 diabetes",    ["pancreas"]),
    ("diabetes mellitus",  "diabetes",           ["pancreas", "adipose tissue", "liver"]),
    ("diabetes",           "diabetes",           ["pancreas", "adipose tissue", "liver"]),
    ("breast cancer",      "breast cancer",      ["breast"]),
    ("lung cancer",        "lung cancer",        ["lung"]),
    ("liver cancer",       "liver cancer",       ["liver"]),
    ("colorectal cancer",  "colorectal cancer",  ["colon"]),
    ("prostate cancer",    "prostate cancer",    ["prostate gland"]),
    ("glioma",             "glioma",             ["brain"]),
    ("glioblastoma",       "glioblastoma",       ["brain"]),
    ("melanoma",           "melanoma",           ["skin"]),
    ("leukemia",           "leukemia",           ["bone marrow", "blood"]),
    ("lymphoma",           "lymphoma",           ["bone marrow"]),
    ("carcinoma",          "carcinoma",          []),
    ("tumor",              "tumor",              []),
    ("tumour",             "tumor",              []),
    ("cancer",             "cancer",             []),
    ("neoplasm",           "cancer",             []),
    ("malignant",          "cancer",             []),
    ("alzheimer",          "Alzheimer's disease", ["brain"]),
    ("parkinson",          "Parkinson's disease", ["brain"]),
    ("asthma",             "asthma",             ["lung"]),
    ("hypertension",       "hypertension",       ["heart", "blood"]),
    ("atherosclerosis",    "atherosclerosis",    ["blood"]),
    ("hepatitis",          "hepatitis",          ["liver"]),
    ("cirrhosis",          "cirrhosis",          ["liver"]),
    ("nephritis",          "nephritis",          ["kidney"]),
    ("osteoarthritis",     "osteoarthritis",     []),
    ("arthritis",          "arthritis",          []),
    ("infection",          "infection",           []),
    ("viral",              "viral infection",     []),
    ("inflammation",       "inflammation",        []),
    ("autoimmune",         "autoimmune disease",  []),
    ("obesity",            "obesity",            ["adipose tissue"]),
    ("fibrosis",           "fibrosis",           []),
    ("ischemia",           "ischemia",           ["heart"]),
    ("muscular dystrophy", "muscular dystrophy", ["skeletal muscle tissue"]),
]

# Legacy alias used by disease→anatomy inference
DISEASE_ANATOMY: dict[str, list[str]] = {
    "type 2 diabetes": ["pancreas", "adipose tissue", "liver"],
    "type 1 diabetes": ["pancreas"],
    "diabetes": ["pancreas", "adipose tissue", "liver"],
    "breast cancer": ["breast"], "lung cancer": ["lung"],
    "liver cancer": ["liver"], "colorectal cancer": ["colon"],
    "prostate cancer": ["prostate gland"],
    "leukemia": ["bone marrow", "blood"],
    "alzheimer": ["brain"], "parkinson": ["brain"],
    "asthma": ["lung"], "hypertension": ["heart", "blood"],
    "hepatitis": ["liver"], "cirrhosis": ["liver"], "nephritis": ["kidney"],
}

_BTO_CELL_HINTS = frozenset(
    {"cell", "cyte", "blast", "neuron", "macrophage", "lymphocyte",
     "fibroblast", "platelet", "monocyte"})

_TEMPORAL_RULES: list[tuple[str, list[str]]] = [
    ("<1s",         ["millisecond", "action potential", "channel gating"]),
    ("1s<1min",     ["signaling", "signalling", "phosphorylation", "binding",
                     "receptor activation"]),
    ("1min<1hr",    ["minute", "metaboli", "enzyme kinetic", "glycoly",
                     "glucose uptake", "insulin secretion"]),
    ("1hr<1day",    ["hour", "circadian", "cell cycle", "gene expression",
                     "transcription"]),
    ("1day<1week",  ["day", "week", "growth", "proliferation",
                     "differentiation", "wound healing"]),
    ("1week<1year", ["month", "year", "chronic", "disease progression",
                     "epidemic", "epidemiolog", "population"]),
]


# -- Helpers -----------------------------------------------------------------

def _clean(text: str) -> str:
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", text)).strip()


def _add(lst: list, item):
    if item and item not in lst:
        lst.append(item)


def _keyword_hits(text: str, table: dict[str, str]) -> set[str]:
    low = text.lower()
    return {v for k, v in table.items() if k in low}


def _dedup(items: list[str]) -> list[str]:
    seen: set[str] = set()
    return [x for x in items if not (x in seen or seen.add(x))]


def _resolve_bto_name(name: str):
    """→ (anat_label, uberon_id, cell_label, cl_id)"""
    if any(h in name.lower() for h in _BTO_CELL_HINTS):
        lab, cid = api.lookup_cl(name)
        return "", "", lab, cid
    lab, uid = api.lookup_uberon(name)
    return lab, uid, "", ""


def _add_anatomy(anat, uberon, name):
    lab, uid = api.lookup_uberon(name)
    _add(anat, lab); _add(uberon, uid)


def _add_cell(cells, cl, name):
    lab, cid = api.lookup_cl(name)
    _add(cells, lab); _add(cl, cid)


def _scan_keywords(text, anat, uberon, cells, cl):
    """Scan text for anatomy/cell keywords, add new hits. Returns set of new terms found."""
    new_anat = _keyword_hits(text, ANATOMY_TERMS) - set(anat)
    new_cells = _keyword_hits(text, CELL_TERMS) - set(cells)
    for name in new_anat:
        _add_anatomy(anat, uberon, name)
    for name in new_cells:
        _add_cell(cells, cl, name)
    return new_anat | new_cells


def _parse_model_annotations(details: dict):
    """Extract structured buckets from modelLevelAnnotations."""
    tissues, go_terms, diseases, organisms, uniprot_ids = [], [], [], [], []
    for a in details.get("modelLevelAnnotations", []):
        q = a.get("qualifier", "")
        res, name, uri, acc = a.get("resource", ""), a.get("name", ""), a.get("uri", ""), a.get("accession", "")
        if q == "bqbiol:occursIn":
            tissues.append(name)
        elif q == "bqbiol:hasTaxon":
            organisms.append(name)
        elif q in ("bqbiol:isVersionOf", "bqbiol:hasPart") and ("Gene Ontology" in res or "go/" in uri.lower()):
            go_terms.append(name)
        elif q == "bqbiol:hasProperty" and ("doid" in uri.lower() or "Human Disease" in res):
            diseases.append(name)
        if "UniProt" in res or "uniprot" in uri.lower():
            _add(uniprot_ids, acc)
    return tissues, go_terms, diseases, organisms, uniprot_ids


# -- Public API --------------------------------------------------------------

def enrich(details: dict, *, use_llm: bool = False) -> dict:
    """Return HRA mapping fields derived from BioModels *details*."""
    tissues, go_terms, diseases, organisms, model_uniprot = _parse_model_annotations(details)
    desc = _clean(details.get("description", ""))
    name = details.get("name", "")
    prov: list[str] = []

    anat, uberon, cells, cl, genes, hgnc = ([] for _ in range(6))
    seen_genes: set[str] = set()

    # 1. SBML file
    sbml = api.fetch_sbml(details)
    all_uniprot = _dedup(model_uniprot + sorted(sbml.uniprot))

    # 1a. Genes from UniProt
    tissue_texts = []
    for uid in all_uniprot[:20]:
        info = api.lookup_uniprot(uid)
        sym, hid = info.get("gene_symbol", ""), info.get("hgnc_id", "")
        if sym and hid and sym not in seen_genes:
            seen_genes.add(sym); genes.append(sym); hgnc.append(hid)
        if info.get("tissue_specificity"):
            tissue_texts.append(info["tissue_specificity"])
    if genes:
        src = "SBML species" if sbml.uniprot else "model-level"
        prov.append(f"Genes: {len(genes)} from {src} UniProt ({', '.join(genes[:5])}{'...' if len(genes)>5 else ''}).")

    # 1b. SBML CL annotations
    for cl_id in sorted(sbml.cl):
        lab, cid = api.ols_lookup(cl_id, "cl", "CL:")
        _add(cells, lab); _add(cl, cid)
    if sbml.cl:
        prov.append(f"Cells: {len(sbml.cl)} from SBML CL annotations.")

    # 1c. SBML BTO annotations
    for bto_id in sorted(sbml.bto):
        lab, _ = api.ols_lookup(bto_id, "bto", "BTO:")
        if lab:
            a, u, c, ci = _resolve_bto_name(lab)
            _add(anat, a); _add(uberon, u); _add(cells, c); _add(cl, ci)
    if sbml.bto:
        prov.append(f"Tissues: {len(sbml.bto)} from SBML BTO.")

    # 2. SBML compartments
    for comp in sbml.compartments:
        lab, uid = api.lookup_uberon(comp)
        if uid:
            _add(anat, lab); _add(uberon, uid)
    if sbml.compartments:
        prov.append(f"Compartments: {', '.join(sbml.compartments)}.")

    # 3. Model-level BTO tissues
    for t in tissues:
        a, u, c, ci = _resolve_bto_name(t)
        _add(anat, a); _add(uberon, u); _add(cells, c); _add(cl, ci)
    if tissues:
        prov.append(f"Model tissues: {', '.join(tissues)}.")

    # 4. PubMed abstract
    pub_link = (details.get("publication") or {}).get("link", "")
    pmid = pub_link.split("pubmed/")[-1] if "pubmed/" in pub_link else ""
    abstract = api.fetch_abstract(pmid)
    if abstract:
        hits = _scan_keywords(abstract, anat, uberon, cells, cl)
        if hits:
            prov.append(f"Abstract keywords: {', '.join(sorted(hits))}.")

    # 5. UniProt tissue specificity (only when anatomy sparse)
    if tissue_texts and len(anat) < 2:
        combined = " ".join(tissue_texts)
        new_a = sorted(_keyword_hits(combined, ANATOMY_TERMS) - set(anat))[:3]
        new_c = sorted(_keyword_hits(combined, CELL_TERMS) - set(cells))[:2]
        for n in new_a:
            _add_anatomy(anat, uberon, n)
        for n in new_c:
            _add_cell(cells, cl, n)
        if new_a or new_c:
            prov.append(f"UniProt tissue: {', '.join((new_a + new_c)[:3])}.")

    # 6. Description + model name keywords
    hits = _scan_keywords(f"{name} {desc}", anat, uberon, cells, cl)
    if hits:
        prov.append(f"Description keywords: {', '.join(sorted(hits))}.")

    # 7. Disease → anatomy inference
    if not anat and diseases:
        for disease in diseases:
            d_low = disease.lower()
            for keyword, _, organs in DISEASE_CONTEXT:
                if keyword in d_low and organs:
                    for org in organs:
                        _add_anatomy(anat, uberon, org)
                    prov.append(f"Anatomy from disease: {disease}.")
                    break

    # Gaps
    if not genes: prov.append("No human genes found.")
    if not anat: prov.append("No anatomy found.")
    if not cells: prov.append("No cell types found.")
    if organisms: prov.append(f"Organism: {', '.join(organisms)}.")
    if diseases: prov.append(f"Disease: {', '.join(diseases)}.")

    # System modeled
    system = "; ".join(go_terms[:3])
    if not system:
        first = desc.split(".")[0].strip()
        if 0 < len(first) < 300:
            system = first

    # Scales
    all_text = f"{name} {desc} {abstract}"
    temporal = _infer_temporal(all_text)
    atlas = _infer_atlas(anat, cells, all_text, bool(sbml.subcellular))

    # Physiological context
    context = _infer_context(diseases, all_text)

    # Confidence
    n = sum(map(bool, [anat, cells, genes, uberon, cl]))
    confidence = "High" if n >= 3 else ("Medium" if n >= 1 else "Low")

    result = {
        "Temporal_Scale": temporal, "Atlas_Scale": atlas,
        "Physiological_Context": context,
        "System_Modeled": system,
        "Anatomical_Structures": ";".join(anat), "Uberon_IDs": ";".join(uberon),
        "Cell_Types": ";".join(cells), "CL_IDs": ";".join(cl),
        "Key_Biomolecules": ";".join(genes), "HGNC_IDs": ";".join(hgnc),
        "Mapping_Confidence": confidence,
        "Mapping_Provenance": " ".join(prov),
    }

    if use_llm:
        from .llm import fill_gaps
        result = fill_gaps(details, result)
    return result


# -- Scale inference ---------------------------------------------------------

def _infer_temporal(text: str) -> str:
    low = text.lower()
    for scale, kws in _TEMPORAL_RULES:
        if any(k in low for k in kws):
            return scale
    return ""


def _infer_context(diseases: list[str], text: str) -> str:
    """Classify model as 'healthy' or a specific disease/condition.

    Checks DOID disease annotations first, then scans text for disease keywords.
    Returns 'healthy' if nothing disease-related is found.
    """
    # Check DOID disease names from annotations
    for disease in diseases:
        d_low = disease.lower()
        for keyword, label, _ in DISEASE_CONTEXT:
            if keyword in d_low:
                return label

    # Scan description/abstract text
    low = text.lower()

    # Explicit healthy markers
    healthy_markers = ["healthy state", "healthy condition", "normal physiology",
                       "healthy volunteer", "healthy subject", "physiological model",
                       "wild type", "wild-type", "normal state"]
    has_healthy = any(m in low for m in healthy_markers)

    # Check disease keywords in text
    for keyword, label, _ in DISEASE_CONTEXT:
        if keyword in low:
            if has_healthy:
                return f"healthy + {label}"
            return label

    return "healthy"


def _infer_atlas(anat, cells, text, has_subcellular) -> str:
    low = text.lower()
    if any(k in low for k in ("whole body", "whole-body", "systemic", "pharmacokinetic", "multi-organ")):
        return "whole body"
    if any(k in low for k in ("population", "epidemic", "epidemiolog")):
        return "population"
    if anat and not cells: return "organ"
    if cells: return "organ" if anat else "cell"
    if has_subcellular: return "cell"
    if any(k in low for k in ("pathway", "signaling", "signalling", "intracellular", "receptor")):
        return "cell"
    if any(k in low for k in ("organ", "tissue")): return "organ"
    return ""
