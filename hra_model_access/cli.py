"""CLI entry point for HRA mapping generation."""

from __future__ import annotations

import argparse
import json
import sys

from . import api
from .enrichment import enrich

FIELDS = [
    "Model_ID", "Model_Name", "Paper_Title", "Paper_Publication_Date",
    "Paper_Link", "DOI", "PubMed_ID", "BioModels_Status", "Model_Type",
    "Temporal_Scale", "Atlas_Scale", "Physiological_Context", "System_Modeled",
    "Anatomical_Structures", "Uberon_IDs", "Cell_Types", "CL_IDs",
    "Key_Biomolecules", "HGNC_IDs", "Mapping_Confidence",
    "Mapping_Provenance", "Primary_Source_URL",
]


def extract_api_fields(model_id: str, details: dict) -> dict:
    """Fields available directly from BioModels metadata."""
    pub = details.get("publication") or {}
    link = pub.get("link", "")

    pubmed_id = link.split("pubmed/")[-1] if "pubmed/" in link else ""
    doi = link.split("doi.org/")[-1] if "doi.org/" in link else ""
    paper_link = (
        f"https://doi.org/{doi}" if doi
        else f"https://pubmed.ncbi.nlm.nih.gov/{pubmed_id}" if pubmed_id
        else link
    )

    year = pub.get("year", "")
    month = pub.get("month", "")
    pub_date = f"{year}-{month:02d}" if isinstance(month, int) and year else str(year)

    approach = details.get("modellingApproach", {}).get("name", "")
    curation = {"CURATED": "Yes - curated", "NON_CURATED": "Non-curated / model entry"}

    return {
        "Model_ID": model_id,
        "Model_Name": details.get("name", ""),
        "Paper_Title": pub.get("title", ""),
        "Paper_Publication_Date": pub_date,
        "Paper_Link": paper_link,
        "DOI": doi,
        "PubMed_ID": pubmed_id,
        "BioModels_Status": curation.get(details.get("curationStatus", ""),
                                         details.get("curationStatus", "")),
        "Model_Type": f"Mathematical Modelling Ontology: {approach}" if approach else "",
        "Primary_Source_URL": f"https://www.biomodels.org/{model_id}",
    }


def process_model(model_id: str, *, use_llm: bool = False) -> dict:
    print(f"  {model_id}...", end=" ", flush=True)
    details = api.fetch_model(model_id)
    row = {f: "" for f in FIELDS}
    row.update(extract_api_fields(model_id, details))
    row.update({k: v for k, v in enrich(details, use_llm=use_llm).items() if v})
    print("done")
    return row


def main(argv: list[str] | None = None):
    p = argparse.ArgumentParser(description="Generate HRA mapping JSON from BioModels")
    p.add_argument("--query", help="BioModels search query (e.g. 'glucose')")
    p.add_argument("--ids", nargs="+", help="Specific BioModels IDs")
    p.add_argument("--max-results", type=int, default=10)
    p.add_argument("--output", default="hra_mapping.json")
    p.add_argument("--no-llm", action="store_true",
                   help="Skip LLM inference, use only API + ontology lookups")
    args = p.parse_args(argv)

    if not args.query and not args.ids:
        p.error("Provide --query or --ids")

    if args.ids:
        model_ids = args.ids
    else:
        print(f"Searching BioModels for '{args.query}'...")
        model_ids = api.search_biomodels(args.query, max_results=args.max_results)
        print(f"Found {len(model_ids)} models")

    if not model_ids:
        print("No models found.")
        sys.exit(1)

    rows = []
    for mid in model_ids:
        try:
            rows.append(process_model(mid, use_llm=not args.no_llm))
        except Exception as e:
            print(f"  ERROR {mid}: {e}")

    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(rows, f, indent=2, ensure_ascii=False)

    # Coverage summary
    print(f"\nWrote {len(rows)} models → {args.output}")
    print("Coverage:")
    for field in FIELDS:
        n = sum(1 for r in rows if r.get(field))
        print(f"  {field}: {n}/{len(rows)} ({n*100//len(rows) if rows else 0}%)")


if __name__ == "__main__":
    main()
