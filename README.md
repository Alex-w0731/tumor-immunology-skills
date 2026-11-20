<p align="center">
  <img src="docs/assets/banner.svg" width="100%" alt="tumor-immunology-skills: evidence-linked agent skills for tumor immunology research. Cited sources, standard-library scripts, unit tests, evals." />
</p>

<p align="center">
  <b>Agent skills that make coding agents careful with tumor-immune data.</b><br />
  Decision rules with checked sources · standard-library scripts with unit tests · evals for every skill
</p>

<p align="center">
  <a href="https://github.com/Alex-w0731/tumor-immunology-skills/actions/workflows/ci.yml"><img src="https://github.com/Alex-w0731/tumor-immunology-skills/actions/workflows/ci.yml/badge.svg" alt="CI" /></a>
  <a href="https://github.com/Alex-w0731/tumor-immunology-skills/actions/workflows/citations.yml"><img src="https://github.com/Alex-w0731/tumor-immunology-skills/actions/workflows/citations.yml/badge.svg" alt="Citations resolve" /></a>
  <a href="https://agentskills.io/specification"><img src="https://img.shields.io/badge/format-Agent%20Skills-7c5cff" alt="Agent Skills format" /></a>
  <img src="https://img.shields.io/badge/scripts-Python%203.9%2B%20stdlib%20only-3776AB?logo=python&logoColor=white" alt="Scripts: Python 3.9+, standard library only" />
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-2E7EA6" alt="MIT license" /></a>
</p>

<p align="center">
  <a href="README.zh-CN.md">中文说明</a> ·
  <a href="#skills">Skills</a> ·
  <a href="#install">Install</a> ·
  <a href="docs/verification.md">How it is verified</a> ·
  <a href="docs/authoring.md">Authoring guide</a>
</p>

---

## Why this exists

A capable coding agent already knows what a Cox model, a UMAP or a deconvolution method is. What goes wrong in tumor immunology is the judgment around them. Left alone, an agent will happily:

- treat 30,000 cells from six patients as thirty thousand independent observations;
- pick the survival cut-point that gives the smallest p-value;
- compare repertoire richness between samples sequenced at different depths;
- read an enrichment score as a cell fraction;
- attach a real citation to a threshold the paper never reported.

Each skill in this repository carries the rule, the reason for it, a tested script for the calculation, and the primary source, so that both the agent and you can check the work. They are written for research use and stay with analysis, statistics, interpretation and reporting.

<a name="skills"></a>

## Skills

<!-- catalog:start -->

**Profile the tumor immune microenvironment**

| Skill | What it makes the agent do | Bundled scripts |
|---|---|---|
| [`tme-bulk-profiling`](skills/tme-bulk-profiling/SKILL.md) | Checks a bulk tumor expression matrix, scores immune signatures and reports immune phenotypes with their limits. | `check_expression_matrix.py` `signature_coverage.py` |
| [`immune-deconvolution-review`](skills/immune-deconvolution-review/SKILL.md) | Chooses, runs and sanity-checks immune deconvolution of bulk tumor transcriptomes, with an NNLS script and a mixture benchmark. | `mixture_benchmark.py` `nnls_deconvolve.py` |

**Single-cell, repertoire and spatial statistics**

| Skill | What it makes the agent do | Bundled scripts |
|---|---|---|
| [`single-cell-patient-level-stats`](skills/single-cell-patient-level-stats/SKILL.md) | Compares groups in single-cell data at the patient level: pseudobulk DE, composition tests, paired designs and few-patient limits. | `composition_test.py` `pseudobulk_aggregate.py` |
| [`tcr-repertoire-metrics`](skills/tcr-repertoire-metrics/SKILL.md) | Computes TCR repertoire diversity, clonality and overlap from clonotype counts, compared at equal sequencing depth. | `repertoire_metrics.py` `repertoire_overlap.py` |

**Outcome and response**

| Skill | What it makes the agent do | Bundled scripts |
|---|---|---|
| [`immune-survival-analysis`](skills/immune-survival-analysis/SKILL.md) | Tests immune features against survival with continuous Cox models, sound endpoints, FDR and no cut-point fishing. | `cutpoint_inflation.py` `fdr_adjust.py` `km_logrank.py` |
| [`response-biomarker-evaluation`](skills/response-biomarker-evaluation/SKILL.md) | Tests whether a biomarker separates immunotherapy responders from non-responders, with honest uncertainty and leakage checks. | `roc_auc.py` `threshold_metrics.py` |

<sub>6 skills, 13 scripts.</sub>

<!-- catalog:end -->

Every skill folder has the same shape:

```text
skills/<name>/
  SKILL.md              what to do, in what order, and why; defaults with source keys; pitfalls; report template
  references/
    sources.md          the bibliography: every key used in the skill, with DOI and PMID
    <topic>.md          detail that is loaded only when needed
  scripts/<tool>.py     Python standard library only, deterministic, documented formulas
  evals/<case>/         prompt, expected outcome and graders (the `claude plugin eval` layout)
tests/<name>/           unit tests for the scripts, with hand-derived expected values
```

<a name="install"></a>

## Install

**Claude Code, as a plugin**

```text
/plugin marketplace add Alex-w0731/tumor-immunology-skills
/plugin install tumor-immunology@tumor-immunology-skills
```

Skills load on their own when a request matches their description. To call one explicitly, type `/tumor-immunology:<skill-name>`. To update later: `/plugin marketplace update tumor-immunology-skills`.

**One skill, by hand**

```bash
git clone https://github.com/Alex-w0731/tumor-immunology-skills.git
mkdir -p ~/.claude/skills && cp -R tumor-immunology-skills/skills/immune-survival-analysis ~/.claude/skills/
```

Each skill is self-contained, so copying its folder is enough.

**Other agents**

The folders follow the open [Agent Skills](https://agentskills.io) format and use only its standard frontmatter fields (`name`, `description`, `license`, `compatibility`, `metadata`). Agents that implement the format can load them: copy a skill folder into the skills directory your agent reads. See the documentation for [Codex](https://developers.openai.com/codex/skills), [Cursor](https://cursor.com/docs/context/skills), [GitHub Copilot](https://docs.github.com/en/copilot/concepts/agents/about-agent-skills), [Gemini CLI](https://geminicli.com/docs/cli/skills/) and [OpenCode](https://opencode.ai/docs/skills).

**Run a script without any agent**

The scripts are ordinary command-line tools with no dependencies:

```bash
python3 skills/immune-survival-analysis/scripts/km_logrank.py --help
```

## How it is verified

| Layer | What is checked | Where |
|---|---|---|
| Sources | Every numeric threshold, definition and attributed claim carries a key that resolves in the skill's `references/sources.md`. A checker resolves every DOI and PMID at Crossref and PubMed and compares year, first author and title with the entry. | [`tools/check_citations.py`](tools/check_citations.py), run by a [workflow](.github/workflows/citations.yml) whenever a bibliography changes |
| Scripts | Standard library only, no network unless declared, `--help` runs, unit tests with expected values derived by hand or from published worked examples. | [`tools/run_tests.py`](tools/run_tests.py), CI on Python 3.9 to 3.13 (Linux, macOS, Windows) |
| Structure | Frontmatter, line budget, link targets, orphan files, eval layout, no machine-specific paths. | [`tools/validate_skills.py`](tools/validate_skills.py) |
| Review | Each skill was drafted by one AI agent and then audited by a second, independent one that re-read the cited sources, recomputed script outputs and repaired what it found. | [`docs/verification.md`](docs/verification.md) |
| Behaviour | Three eval cases per skill with objective graders, meant to be run with and without the skill. | `skills/<name>/evals/` |

The machine checks catch references that do not exist or that point at a different paper. They cannot prove that a paper supports the sentence it is cited for; that depends on review, and review can miss things. [`docs/verification.md`](docs/verification.md) says plainly what has and has not been checked. If you find an error, please open an issue with the source: corrections are the most valuable contribution here.

## Scope and limits

- **Research use only.** Nothing here is clinical decision support.
- **Analysis side only.** Data analysis, statistics, interpretation and reporting. No wet-lab protocols.
- **Defaults, not laws.** Thresholds are given with their provenance and labelled as published, conventional or heuristic. Date-sensitive facts carry an "as of" stamp.
- **Synthetic examples are labelled synthetic** and are never presented as findings.

## Related

- [`tmelens`](https://github.com/Alex-w0731/tmelens): a transparent Python toolkit to profile the tumour immune microenvironment from bulk transcriptomes. Several skills here use it when it is installed and fall back to their own scripts when it is not.
- The [Agent Skills specification](https://agentskills.io/specification).

## Contributing

Corrections and new skills are welcome. Read [CONTRIBUTING.md](CONTRIBUTING.md) and the [authoring guide](docs/authoring.md) first: the bar is that a skill changes what an agent does, cites what it claims, and tests what it computes.

```bash
python3 tools/validate_skills.py     # structure, evidence keys, scripts, evals
python3 tools/run_tests.py           # unit tests of every bundled script
python3 tools/build_catalog.py       # refresh the catalog tables above
```

## Citation

If these skills contribute to your work, cite the repository (see [`CITATION.cff`](CITATION.cff)) together with the primary papers listed in the relevant skill's `references/sources.md`.

## License

[MIT](LICENSE) © Alex Wang
