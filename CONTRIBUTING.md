# Contributing

Thank you for helping. The most valuable contribution is a correction: a claim that the cited
source does not support, a threshold that has been superseded, a script that computes the wrong
thing. Please open an issue and include the source.

## The bar for a skill

A skill belongs here when all three hold:

1. **It changes what an agent does.** It carries decision rules, pitfalls, defaults with
   provenance or a tested calculation that a capable general agent would otherwise get wrong,
   skip, or reinvent badly. Restating a library's documentation does not qualify.
2. **It cites what it claims.** Every numeric threshold, definition, gene list or attributed
   statement has a source key that resolves in `references/sources.md`, and the source says it.
3. **It tests what it computes.** Scripts use the Python standard library only and have unit
   tests with expected values derived by hand or from a published worked example.

The full rules are in the [authoring guide](docs/authoring.md). Scope is the analysis side of
tumor immunology research: data analysis, statistics, interpretation and reporting. Wet-lab
protocols and anything that reads as clinical advice are out of scope.

## Checks to run before opening a pull request

```bash
python3 tools/validate_skills.py          # structure, evidence keys, scripts, evals
python3 tools/run_tests.py                # unit tests of every bundled script
python3 tools/build_catalog.py            # refresh the catalog tables in both READMEs
python3 tools/check_citations.py          # resolve every DOI and PMID (needs network)
```

No dependencies need to be installed. CI runs the first three on Python 3.10 to 3.14; a separate
workflow runs the citation check whenever a bibliography changes.

## Review checklist

A reviewer (human or agent) should be able to answer yes to each of these:

- Is every cited source real, and does it support the sentence it is attached to? (Read the
  abstract, and the full text when a specific number is attributed.)
- Are published thresholds, common conventions and the skill's own heuristics told apart?
- Do date-sensitive statements carry an "as of" stamp?
- Were at least two outputs of each script recomputed independently?
- Would an agent that received only `SKILL.md` know what to do at each decision point?
- Would a capable agent without the skill plausibly fail the eval graders?

## Versioning

Each skill has its own `metadata.version`; bump it when the skill changes. The plugin version in
`.claude-plugin/plugin.json` and `CHANGELOG.md` track releases of the collection.

## Conduct

Be precise and be kind. Disagreements about evidence are settled by reading the source together.
