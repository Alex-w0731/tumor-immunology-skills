#!/usr/bin/env python3
"""Regenerate the skill catalog tables in README.md and README.zh-CN.md.

The tables are built from each skill's frontmatter (`metadata.summary`, `metadata.summary-zh`)
and its scripts/ folder, and written between the markers

    <!-- catalog:start -->  ...  <!-- catalog:end -->

so the READMEs can never drift from the skills that are actually in the repository.

usage:
  python3 tools/build_catalog.py           # rewrite both READMEs
  python3 tools/build_catalog.py --check   # exit 1 if a README is out of date (used in CI)
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
START, END = "<!-- catalog:start -->", "<!-- catalog:end -->"

# Display order. A skill that is not listed here is appended under "More".
GROUPS = [
    ("Profile the tumor immune microenvironment", "刻画肿瘤免疫微环境",
     ["tme-bulk-profiling", "immune-deconvolution-review", "immune-signature-curation"]),
    ("Single-cell, repertoire and spatial statistics", "单细胞、受体库与空间统计",
     ["single-cell-patient-level-stats", "tcr-repertoire-metrics", "spatial-neighborhood-statistics"]),
    ("Outcome and response", "生存结局与疗效",
     ["immune-survival-analysis", "response-biomarker-evaluation"]),
    ("Evidence", "证据核查",
     ["citation-evidence-check"]),
]
HEADERS = {
    "en": ("Skill", "What it makes the agent do", "Bundled scripts"),
    "zh": ("技能", "它让智能体做对什么", "附带脚本"),
}


def _validator():
    spec = importlib.util.spec_from_file_location("validate_skills", ROOT / "tools" / "validate_skills.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["validate_skills"] = module
    spec.loader.exec_module(module)
    return module


def load_skills() -> dict[str, dict]:
    v = _validator()
    skills: dict[str, dict] = {}
    for folder in sorted((ROOT / "skills").iterdir()):
        md = folder / "SKILL.md"
        if not md.is_file():
            continue
        front, _ = v.split_frontmatter(md.read_text(encoding="utf-8"))
        meta = v.parse_simple_yaml(front or "")
        extra = meta.get("metadata") if isinstance(meta.get("metadata"), dict) else {}
        scripts = sorted(p.name for p in (folder / "scripts").glob("*.py")) if (folder / "scripts").is_dir() else []
        skills[folder.name] = {
            "summary": str(extra.get("summary", "")).strip(),
            "summary-zh": str(extra.get("summary-zh", "")).strip(),
            "scripts": scripts,
        }
    return skills


def cell(text: str) -> str:
    return text.replace("|", "\\|").replace("\n", " ")


def render(skills: dict[str, dict], lang: str) -> str:
    head = HEADERS[lang]
    key = "summary" if lang == "en" else "summary-zh"
    listed = {name for _, _, names in GROUPS for name in names}
    groups = [(en if lang == "en" else zh, names) for en, zh, names in GROUPS]
    rest = [n for n in skills if n not in listed]
    if rest:
        groups.append(("More" if lang == "en" else "更多", rest))
    out: list[str] = []
    for title, names in groups:
        present = [n for n in names if n in skills]
        if not present:
            continue
        out.append(f"**{title}**")
        out.append("")
        out.append(f"| {head[0]} | {head[1]} | {head[2]} |")
        out.append("|---|---|---|")
        for name in present:
            info = skills[name]
            scripts = " ".join(f"`{s}`" for s in info["scripts"]) or "-"
            out.append(f"| [`{name}`](skills/{name}/SKILL.md) | {cell(info[key] or info['summary'])} | {scripts} |")
        out.append("")
    count = len(skills)
    scripts_total = sum(len(s["scripts"]) for s in skills.values())
    footer = (f"{count} skills, {scripts_total} scripts." if lang == "en" else f"共 {count} 个技能、{scripts_total} 个脚本。")
    out.append(f"<sub>{footer}</sub>")
    return "\n".join(out)


def splice(text: str, block: str) -> str:
    if START not in text or END not in text:
        raise SystemExit("catalog markers not found")
    before, rest = text.split(START, 1)
    _, after = rest.split(END, 1)
    return f"{before}{START}\n\n{block}\n\n{END}{after}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build the README skill catalogs.")
    parser.add_argument("--check", action="store_true", help="verify that the READMEs are up to date")
    args = parser.parse_args(argv)
    skills = load_skills()
    stale = []
    for filename, lang in (("README.md", "en"), ("README.zh-CN.md", "zh")):
        path = ROOT / filename
        if not path.is_file():
            if args.check:
                stale.append(f"{filename} (missing)")
            continue
        current = path.read_text(encoding="utf-8")
        updated = splice(current, render(skills, lang))
        if updated != current:
            if args.check:
                stale.append(filename)
            else:
                path.write_text(updated, encoding="utf-8")
                print(f"updated {filename}")
    if args.check:
        if stale:
            print("catalog out of date (run python3 tools/build_catalog.py):", ", ".join(stale))
            return 1
        print(f"catalogs up to date ({len(skills)} skills)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
