"""CLI for evaluating HRA mapping JSONs with COPASI simulation."""

from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None):
    p = argparse.ArgumentParser(
        description="Test which BioModels in a mapping JSON are runnable via COPASI")
    p.add_argument("input", help="Input mapping JSON (from generate-hra-mapping)")
    p.add_argument("-o", "--output", help="Output JSON with Runs_In_COPASI column "
                   "(default: overwrites input)")
    p.add_argument("--html", help="Write HTML report with time-course plots")
    p.add_argument("--max-models", type=int, default=None,
                   help="Only test the first N models")
    args = p.parse_args(argv)

    try:
        from .simulate import evaluate_mapping
    except ImportError as e:
        print(f"Missing dependency: {e}\n"
              f"Install with: pip install -e '.[simulate]'", file=sys.stderr)
        sys.exit(1)

    output = args.output or args.input
    evaluate_mapping(args.input, output, max_models=args.max_models,
                     html_report=args.html)


if __name__ == "__main__":
    main()
