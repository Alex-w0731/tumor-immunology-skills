# How these skills are verified

This page says what has been checked, how, and what has not. It is meant to be read before you
rely on a skill.

## How the skills were made

The skills were drafted by AI coding agents (Claude) working from a written brief and the
[authoring guide](authoring.md), under the direction of the repository owner. Each skill then
went to a second, independent agent whose only job was to find and repair errors: it re-resolved
every cited source, read the abstract (and the open-access full text where a number was
attributed), recomputed script outputs on small inputs without reusing the script's code, and
tightened the instructions and eval graders. Commits carry a co-author trailer for the model.

AI-drafted and AI-reviewed is not the same as reviewed by a domain expert line by line. Treat
the skills as carefully sourced working notes: check the primary source before you build on a
specific number.

## What is checked by machine

| Check | Tool | Guarantees | Does not guarantee |
|---|---|---|---|
| Every source key used in a skill is defined, and every entry has a DOI, PMID or URL | `tools/validate_skills.py` | No dangling citations | That the source is the right one |
| Every DOI and PMID resolves, and year, first author and title match the entry | `tools/check_citations.py` (workflow on every change to a `sources.md`, and on demand) | The reference exists and is the paper described | That the paper supports the sentence it is cited for |
| Scripts import only the standard library, expose `--help`, and have a test file | `tools/validate_skills.py` | Portability and a minimal contract | Correctness |
| Unit tests with hand-derived expected values, edge cases and a command-line smoke test | `tools/run_tests.py`, CI matrix | The implemented formula gives the expected numbers on the tested inputs | Correctness on inputs unlike the tests |
| Layout, line budgets, link targets, eval case format, no machine-specific paths or e-mail addresses | `tools/validate_skills.py` | The files load and are self-contained | Quality of the writing |

## What depends on review

- Whether each claim is supported by its source.
- Whether a default is still current. Date-sensitive statements carry an "as of" stamp; assume
  anything regulatory or guideline-based may have moved since.
- Whether the advice fits your design. The skills give defaults and the reasoning behind them
  so that you can depart from them knowingly.

## Evals

Each skill ships three cases under `skills/<name>/evals/` in the layout that
`claude plugin eval` reads: a realistic prompt, the expected outcome, and three to six objective
graders. They are designed to be run twice, with and without the skill, so that the difference
can be measured rather than asserted. Results, when they have been run, are recorded in
[`evals.md`](evals.md) with the date, the model and the pass rates of both configurations.

## Reporting a problem

Open an issue with: the skill, the sentence or the script output, what you believe is correct,
and the source. Corrections are credited in the changelog.
