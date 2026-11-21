# Changelog

All notable changes to this collection are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the collection uses
[semantic versioning](https://semver.org/).

## [0.1.0] - 2026-11-20

First public release.

### Added

- Nine analysis skills for tumor immunology research, each with a `SKILL.md`, a checked
  bibliography, standard-library scripts with unit tests, and three eval cases:
  `tme-bulk-profiling`, `immune-deconvolution-review`, `immune-signature-curation`,
  `single-cell-patient-level-stats`, `tcr-repertoire-metrics`, `spatial-neighborhood-statistics`,
  `immune-survival-analysis`, `response-biomarker-evaluation` and `citation-evidence-check`.
- Each skill was drafted by one agent and audited by an independent one that re-checked the
  cited sources and recomputed script outputs; the reviews corrected misattributed claims, wrong
  small-sample arithmetic and several script bugs before release.
- First with-skill / without-skill eval run (one case per skill), reported in `docs/evals.md`.
- `tools/validate_skills.py`: structural, evidence-key, script and eval checks.
- `tools/run_tests.py`: runs every skill's unit tests in isolation.
- `tools/check_citations.py`: resolves every cited DOI and PMID at Crossref and PubMed and
  compares the metadata with the entry.
- `tools/build_catalog.py`: keeps the README catalogs in sync with the skills.
- Claude Code plugin and marketplace manifests; CI on Python 3.10 to 3.14 (Linux), plus macOS and Windows.
- Scripts switch standard output to UTF-8, so reports print on Windows consoles with legacy code pages.
