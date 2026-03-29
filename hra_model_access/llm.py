"""Optional LLM-based gap-filling for HRA fields."""

from __future__ import annotations

import json
import re

_FIELD_HELP = {
    "Temporal_Scale": 'time scale (e.g., "1s<1min", "1min<1hr", "1week<1year")',
    "Atlas_Scale": 'spatial scale (e.g., "cell", "organ", "whole body")',
    "Physiological_Context": '"healthy" or disease name (e.g., "tumor", "diabetes", "Parkinson\'s disease")',
    "System_Modeled": "brief description of what the model simulates",
    "Anatomical_Structures": 'semicolon-separated (e.g., "Pancreas;Blood circulatory system")',
    "Uberon_IDs": 'semicolon-separated (e.g., "UBERON:0001264")',
    "Cell_Types": 'semicolon-separated (e.g., "Pancreatic beta cell")',
    "CL_IDs": 'semicolon-separated (e.g., "CL:0000169")',
    "Key_Biomolecules": 'semicolon-separated gene symbols (e.g., "INS;INSR")',
    "HGNC_IDs": 'semicolon-separated (e.g., "HGNC:6081")',
    "Mapping_Confidence": '"High", "Medium", or "Low"',
    "Mapping_Provenance": "one paragraph explaining what sources were used and why gaps remain",
}


def fill_gaps(details: dict, current: dict) -> dict:
    """Call Claude to fill empty fields in *current*; returns updated copy."""
    empty = [f for f in _FIELD_HELP if not current.get(f)]
    if not empty:
        return current

    try:
        import anthropic
    except ImportError:
        print("  WARNING: anthropic not installed — skipping LLM inference")
        return current

    desc = re.sub(r"<[^>]+>", " ", details.get("description", ""))[:3000]
    anns = "\n".join(
        f"  {a.get('qualifier','')}: {a.get('name','')} ({a.get('uri','')})"
        for a in details.get("modelLevelAnnotations", [])
    )
    filled = json.dumps({k: v for k, v in current.items() if v}, indent=2)
    fields = "\n".join(f'- "{f}": {_FIELD_HELP[f]}' for f in empty)

    prompt = (
        f"Analyze this BioModels entry and return a JSON object with ONLY these missing HRA fields:\n"
        f"{fields}\n\n"
        f"Model: {details.get('name','')}\n"
        f"Description: {desc}\n"
        f"Annotations:\n{anns}\n\n"
        f"Already extracted (do not override):\n{filled}\n\n"
        f"Return ONLY valid JSON, no markdown fences."
    )

    try:
        client = anthropic.Anthropic()
        msg = client.messages.create(
            model="claude-sonnet-4-20250514",
            max_tokens=1024,
            messages=[{"role": "user", "content": prompt}],
        )
        text = msg.content[0].text.strip()
        if text.startswith("```"):
            text = text.split("\n", 1)[1].rsplit("```", 1)[0].strip()
        llm = json.loads(text)
    except Exception as e:
        print(f"  WARNING: LLM inference failed: {e}")
        return current

    out = dict(current)
    for k, v in llm.items():
        if k in _FIELD_HELP and not out.get(k):
            out[k] = v
    return out
