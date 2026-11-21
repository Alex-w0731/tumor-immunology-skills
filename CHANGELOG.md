# Changelog

All notable changes to this collection are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the collection uses
[semantic versioning](https://semver.org/).

## [0.1.1] - 2026-11-20

Second iteration, driven by the first blind evaluation.

### Changed

- Eval graders in all 27 cases were audited. Statements the blind grader judged
  non-discriminating (passed for reasons unrelated to the skill) or ambiguous were rewritten or
  removed; statements now say whether a value must appear in the reply or may sit in an output
  file, and every number in a grader was recomputed with the bundled scripts.
- Each skill gained the rules where the evaluation showed a gap, for example: offer the paper
  figure from the valid analysis whatever grouping was requested (`immune-survival-analysis`);
  do not switch tests because an exact floor is above 0.05 (`single-cell-patient-level-stats`);
  a spreadsheet-damage candidate is a lookup, not a repair (`immune-signature-curation`);
  saying a paper does not report something needs its full text (`citation-evidence-check`);
  report observed and expected counts under each null per sample (`spatial-neighborhood-statistics`);
  a benchmark is evidence only for the range it spans (`immune-deconvolution-review`); keep a
  per-patient inclusion table (`response-biomarker-evaluation`); rerun with the other clonotype key
  when the key is a judgement call (`tcr-repertoire-metrics`); score CD8 in tumor and stroma
  separately when testing an exclusion statement (`tme-bulk-profiling`).
- Every skill's `metadata.version` is now 0.1.1. The iteration-1 results in `docs/evals.md` were
  obtained with the 0.1.0 graders and are kept as they were.

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
