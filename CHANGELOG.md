# Changelog

All notable changes to this collection are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the collection uses
[semantic versioning](https://semver.org/).

## [0.1.0] - 2026-11-20

First public release.

### Added

- The first set of analysis skills for tumor immunology research, each with a `SKILL.md`,
  cited references, standard-library scripts, unit tests and three eval cases. See the catalog
  in the README for the current list.
- `tools/validate_skills.py`: structural, evidence-key, script and eval checks.
- `tools/run_tests.py`: runs every skill's unit tests in isolation.
- `tools/check_citations.py`: resolves every cited DOI and PMID at Crossref and PubMed and
  compares the metadata with the entry.
- `tools/build_catalog.py`: keeps the README catalogs in sync with the skills.
- Claude Code plugin and marketplace manifests; CI on Python 3.9 to 3.13.
