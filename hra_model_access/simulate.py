"""Run BioModels SBML models via COPASI (basico) and return time-series results.

Requires: pip install copasi-basico python-libsedml

Usage:
    from hra_model_access.simulate import simulate_model
    result = simulate_model("BIOMD0000000356")
    # result = {"model_id": ..., "time": [...], "species": {"INS": [...], ...}, ...}
"""

from __future__ import annotations

import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree as ET

from . import api


@dataclass
class SimResult:
    model_id: str
    model_name: str
    time: list[float]
    species: dict[str, list[float]]  # species_name → concentration timeseries
    duration: float
    n_points: int
    time_unit: str


def simulate_model(
    model_id: str,
    *,
    duration: float | None = None,
    n_points: int | None = None,
    output_dir: str | None = None,
) -> SimResult:
    """Fetch a BioModels model, run a COPASI time-course, return results.

    If *duration* or *n_points* are None, they are read from the SED-ML file
    bundled with the model (if available), otherwise defaults to 10s / 100 points.
    """
    try:
        from basico import load_model, run_time_course, get_species
    except ImportError:
        raise ImportError("copasi-basico is required: pip install copasi-basico")

    details = api.fetch_model(model_id)
    sbml_path, sed_duration, sed_points = _fetch_sbml(details, output_dir)

    if duration is None:
        duration = sed_duration
    if n_points is None:
        n_points = sed_points

    # Load and run
    dm = load_model(sbml_path)
    if dm is None:
        raise RuntimeError(f"COPASI could not load {sbml_path}")

    tc = _run_with_fallback(dm, duration, max(n_points - 1, 1))

    # Extract species name mapping (SBML ID → display name)
    spec_df = get_species(model=dm)
    sbml_to_name = {}
    if spec_df is not None and "sbml_id" in spec_df.columns:
        sbml_to_name = {row["sbml_id"]: name for name, row in spec_df.iterrows()}

    # Build result
    species = {}
    for col in tc.columns:
        display = sbml_to_name.get(col, col)
        species[display] = tc[col].tolist()

    # Infer time unit from model
    time_unit = _detect_time_unit(sbml_path)

    return SimResult(
        model_id=model_id,
        model_name=details.get("name", model_id),
        time=tc.index.tolist(),
        species=species,
        duration=duration,
        n_points=n_points,
        time_unit=time_unit,
    )


def simulate_from_mapping(mapping_json: str, *, max_models: int | None = None):
    """Run simulations for models in a mapping JSON file. Yields (model_id, SimResult | error)."""
    import json
    with open(mapping_json) as f:
        entries = json.load(f)

    for i, entry in enumerate(entries):
        if max_models and i >= max_models:
            break
        model_id = entry["Model_ID"]
        try:
            result = simulate_model(model_id)
            yield model_id, result
        except Exception as e:
            yield model_id, e


def evaluate_mapping(input_path: str, output_path: str, *,
                     max_models: int | None = None, html_report: str | None = None):
    """Load a mapping JSON, test each model with COPASI, add Runs_In_COPASI column.

    If *html_report* is given, also write an HTML page with time-course plots.
    """
    import json

    with open(input_path) as f:
        entries = json.load(f)

    if max_models:
        entries = entries[:max_models]

    sim_results: list[tuple[dict, SimResult | None]] = []

    for entry in entries:
        model_id = entry["Model_ID"]
        print(f"  {model_id}...", end=" ", flush=True)
        try:
            result = simulate_model(model_id)
            entry["Runs_In_COPASI"] = True
            entry["Simulation_Duration"] = result.duration
            entry["Simulation_Points"] = result.n_points
            entry["Simulation_Species_Count"] = len(result.species)
            sim_results.append((entry, result))
            print(f"OK ({len(result.species)} species, {result.duration}{result.time_unit})")
        except Exception as e:
            entry["Runs_In_COPASI"] = False
            entry["Simulation_Error"] = str(e)
            sim_results.append((entry, None))
            print(f"FAIL: {e}")

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(entries, f, indent=2, ensure_ascii=False)

    n_ok = sum(1 for e in entries if e.get("Runs_In_COPASI"))
    print(f"\nWrote {len(entries)} models → {output_path}")
    print(f"Runnable: {n_ok}/{len(entries)} ({n_ok*100//len(entries) if entries else 0}%)")

    if html_report:
        _write_html_report(sim_results, html_report)
        print(f"Report: {html_report}")


def _write_html_report(sim_results: list[tuple[dict, SimResult | None]], path: str):
    """Generate an HTML report with summary table and time-course plots."""
    n_total = len(sim_results)
    n_ok = sum(1 for _, r in sim_results if r is not None)

    cards = []
    for entry, result in sim_results:
        mid = entry["Model_ID"]
        name = entry.get("Model_Name", mid)
        status = "OK" if result else "FAIL"
        color = "#2d7d46" if result else "#c0392b"

        if result:
            svg = _make_svg_plot(result)
            info = (f"<b>Duration:</b> {result.duration} {result.time_unit} · "
                    f"<b>Species:</b> {len(result.species)} · "
                    f"<b>Points:</b> {result.n_points}")
        else:
            svg = f'<div class="fail-msg">{entry.get("Simulation_Error", "Unknown error")}</div>'
            info = ""

        context = entry.get("Physiological_Context", "")
        anatomy = entry.get("Anatomical_Structures", "")
        genes = entry.get("Key_Biomolecules", "")

        cards.append(f"""
        <div class="card">
          <div class="card-header">
            <span class="status" style="background:{color}">{status}</span>
            <a href="https://www.biomodels.org/{mid}" target="_blank"><b>{mid}</b></a>
            — {_esc(name)}
          </div>
          <div class="card-meta">
            {info}
            {f'<br><b>Context:</b> {_esc(context)}' if context else ''}
            {f' · <b>Anatomy:</b> {_esc(anatomy)}' if anatomy else ''}
            {f' · <b>Genes:</b> {_esc(genes)}' if genes else ''}
          </div>
          <div class="card-plot">{svg}</div>
        </div>""")

    html = f"""<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>HRA Model Evaluation Report</title>
<style>
  body {{ font-family: -apple-system, system-ui, sans-serif; max-width: 1200px;
         margin: 0 auto; padding: 20px; background: #f5f5f5; }}
  h1 {{ margin-bottom: 5px; }}
  .summary {{ background: #fff; padding: 15px; border-radius: 8px; margin-bottom: 20px;
              box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
  .card {{ background: #fff; border-radius: 8px; margin-bottom: 12px; padding: 16px;
           box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}
  .card-header {{ font-size: 14px; margin-bottom: 6px; }}
  .card-meta {{ font-size: 12px; color: #666; margin-bottom: 10px; }}
  .card-plot {{ overflow-x: auto; }}
  .status {{ display: inline-block; color: #fff; padding: 2px 8px; border-radius: 4px;
             font-size: 11px; font-weight: bold; margin-right: 6px; }}
  .fail-msg {{ color: #c0392b; font-size: 13px; padding: 10px 0; }}
  a {{ color: #2563eb; text-decoration: none; }}
  a:hover {{ text-decoration: underline; }}
  svg {{ display: block; }}
</style>
</head><body>
<h1>HRA Model Evaluation Report</h1>
<div class="summary">
  <b>{n_ok}/{n_total}</b> models ran successfully in COPASI
  ({n_ok*100//n_total if n_total else 0}%)
</div>
{''.join(cards)}
</body></html>"""

    Path(path).write_text(html, encoding="utf-8")


def _make_svg_plot(result: SimResult, width: int = 700, height: int = 200) -> str:
    """Generate a simple inline SVG time-course plot."""
    if not result.time or not result.species:
        return ""

    t = result.time
    t_min, t_max = min(t), max(t)
    if t_max == t_min:
        return ""

    margin_l, margin_r, margin_t, margin_b = 50, 20, 10, 30
    pw = width - margin_l - margin_r
    ph = height - margin_t - margin_b

    # Find global y range
    all_vals = [v for series in result.species.values() for v in series if v == v]  # skip NaN
    if not all_vals:
        return ""
    y_min = min(0, min(all_vals))
    y_max = max(all_vals) * 1.1 or 1.0

    def tx(v):
        return margin_l + (v - t_min) / (t_max - t_min) * pw

    def ty(v):
        return margin_t + ph - (v - y_min) / (y_max - y_min) * ph

    # Colors
    palette = ["#2563eb", "#dc2626", "#16a34a", "#d97706", "#7c3aed",
               "#0891b2", "#be185d", "#65a30d", "#a855f7", "#f59e0b"]

    lines = []
    legend = []
    for i, (name, vals) in enumerate(list(result.species.items())[:10]):
        color = palette[i % len(palette)]
        points = []
        for j, v in enumerate(vals):
            if v == v:  # skip NaN
                points.append(f"{tx(t[j]):.1f},{ty(v):.1f}")
        if points:
            lines.append(f'<polyline points="{" ".join(points)}" '
                         f'fill="none" stroke="{color}" stroke-width="1.5" opacity="0.8"/>')
            legend.append(f'<text x="{width - margin_r + 5}" y="{margin_t + 12 + i * 14}" '
                          f'font-size="10" fill="{color}">{_esc(name[:20])}</text>')

    # Axes
    axes = (
        f'<line x1="{margin_l}" y1="{margin_t + ph}" x2="{margin_l + pw}" '
        f'y2="{margin_t + ph}" stroke="#ccc" stroke-width="1"/>'
        f'<line x1="{margin_l}" y1="{margin_t}" x2="{margin_l}" '
        f'y2="{margin_t + ph}" stroke="#ccc" stroke-width="1"/>'
        f'<text x="{margin_l - 5}" y="{margin_t + 10}" font-size="9" '
        f'text-anchor="end" fill="#999">{y_max:.2g}</text>'
        f'<text x="{margin_l - 5}" y="{margin_t + ph}" font-size="9" '
        f'text-anchor="end" fill="#999">{y_min:.2g}</text>'
        f'<text x="{margin_l}" y="{margin_t + ph + 15}" font-size="9" fill="#999">{t_min:.2g}</text>'
        f'<text x="{margin_l + pw}" y="{margin_t + ph + 15}" font-size="9" '
        f'text-anchor="end" fill="#999">{t_max:.2g} {result.time_unit}</text>'
    )

    legend_width = 150 if legend else 0
    total_w = width + legend_width

    return (f'<svg width="{total_w}" height="{height}" xmlns="http://www.w3.org/2000/svg">'
            f'{axes}{"".join(lines)}{"".join(legend)}</svg>')


def _esc(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _fetch_sbml(details: dict, output_dir: str | None) -> tuple[str, float, int]:
    """Download SBML (and optionally SED-ML) from BioModels, return (path, duration, n_points)."""
    model_id = details.get("publicationId") or details.get("submissionId", "")
    all_files = (details.get("files", {}).get("main", [])
                 + details.get("files", {}).get("additional", []))

    sbml_files = [f for f in all_files if f["name"].endswith((".xml", ".sbml"))]
    sedml_files = [f for f in all_files if f["name"].endswith(".sedml")]

    if not sbml_files:
        raise ValueError(f"{model_id}: no SBML file found in BioModels entry")

    # Download SBML
    dest = output_dir or tempfile.mkdtemp(prefix=f"hra_sim_{model_id}_")
    os.makedirs(dest, exist_ok=True)

    sbml_content = api._get(
        f"{api.BIOMODELS}/model/download/{model_id}",
        params={"filename": sbml_files[0]["name"]}, mode="bytes")
    if not sbml_content:
        raise ValueError(f"{model_id}: failed to download SBML")

    sbml_path = os.path.join(dest, sbml_files[0]["name"])
    Path(sbml_path).write_bytes(sbml_content)

    # Try to get duration/n_points from SED-ML
    duration, n_points = 10.0, 100
    if sedml_files:
        sedml_content = api._get(
            f"{api.BIOMODELS}/model/download/{model_id}",
            params={"filename": sedml_files[0]["name"]}, mode="bytes")
        if sedml_content:
            sedml_path = os.path.join(dest, sedml_files[0]["name"])
            Path(sedml_path).write_bytes(sedml_content)
            d, n = _parse_sedml(sedml_path)
            if d and n:
                duration, n_points = d, n

    return sbml_path, duration, n_points


def _parse_sedml(path: str) -> tuple[float | None, int | None]:
    """Extract (duration, n_points) from first UniformTimeCourse in SED-ML."""
    try:
        import libsedml
    except ImportError:
        return None, None
    doc = libsedml.readSedMLFromFile(path)
    if doc is None:
        return None, None
    for i in range(doc.getNumSimulations()):
        sim = doc.getSimulation(i)
        if hasattr(sim, "getOutputEndTime") and hasattr(sim, "getNumberOfPoints"):
            duration = float(sim.getOutputEndTime()) - float(sim.getOutputStartTime())
            return duration, int(sim.getNumberOfPoints())
    return None, None


def _run_with_fallback(dm, duration: float, intervals: int):
    """Run COPASI time-course with automatic solver fallback on NaN."""
    from basico import run_time_course

    kwargs = dict(start_time=0.0, duration=duration, intervals=intervals,
                  update_model=True, use_sbml_id=True, model=dm)

    tc = run_time_course(**kwargs)

    if tc.isnull().any().any():
        tc = run_time_course(**{**kwargs, "max_steps": 500000})

    if tc.isnull().any().any():
        tc = run_time_course(**{**kwargs, "method": "radau5"})

    if tc.isnull().any().any():
        tc = run_time_course(**{**kwargs, "method": "radau5",
                                "r_tol": 1e-3, "a_tol": 1e-6, "max_steps": 10000000})
    return tc


def _detect_time_unit(sbml_path: str) -> str:
    """Detect time unit from SBML model.

    Checks (in order):
    1. model timeUnits attribute
    2. unitDefinition with id/name containing 'time'
    3. Falls back to 's'
    """
    try:
        tree = ET.parse(sbml_path)
        root = tree.getroot()
    except Exception:
        return "s"

    model_time_units = ""
    unit_defs: dict[str, float] = {}  # id → seconds multiplier

    for elem in root.iter():
        tag = elem.tag.rsplit("}", 1)[-1]

        if tag == "model":
            model_time_units = elem.get("timeUnits", "")

        elif tag == "unitDefinition":
            uid = elem.get("id", "")
            # Parse the contained <unit> elements to compute seconds multiplier
            total_seconds = None
            for child in elem.iter():
                ctag = child.tag.rsplit("}", 1)[-1]
                if ctag == "unit" and child.get("kind", "") == "second":
                    mult = float(child.get("multiplier", "1"))
                    scale = int(child.get("scale", "0"))
                    exp = float(child.get("exponent", "1"))
                    total_seconds = mult * (10 ** scale) ** exp
            if total_seconds is not None:
                unit_defs[uid] = total_seconds

    # 1. Check explicit timeUnits attribute
    if model_time_units:
        label = _seconds_to_label(unit_defs.get(model_time_units))
        if label:
            return label
        # Direct name match
        for name, label in [("second", "s"), ("s", "s"), ("minute", "min"),
                            ("min", "min"), ("hour", "h"), ("h", "h"),
                            ("day", "day"), ("d", "day")]:
            if model_time_units == name:
                return label

    # 2. Check unitDefinition named 'time'
    for uid, seconds in unit_defs.items():
        if "time" in uid.lower():
            label = _seconds_to_label(seconds)
            if label:
                return label

    # 3. If there's only one unit def with seconds, use it
    if len(unit_defs) == 1:
        label = _seconds_to_label(next(iter(unit_defs.values())))
        if label:
            return label

    return "s"


def _seconds_to_label(seconds: float | None) -> str:
    """Convert a seconds multiplier to a human label."""
    if seconds is None:
        return ""
    if abs(seconds - 1.0) < 0.1:
        return "s"
    if abs(seconds - 60.0) < 1.0:
        return "min"
    if abs(seconds - 3600.0) < 60.0:
        return "h"
    if abs(seconds - 86400.0) < 600.0:
        return "day"
    # Uncommon but handle weeks/years
    if abs(seconds - 604800.0) < 3600.0:
        return "week"
    return ""
