# Authoring guide

Read this fully before writing or reviewing a skill.

## Who reads a skill

An AI coding agent (Claude Code, Codex, and similar) that is helping a tumour-immunology
researcher. The agent is already a strong generalist that knows what scRNA-seq or a Cox model
is. A skill earns its place only by changing what the agent does: it carries the decision
rules, thresholds with provenance, pitfalls, reporting standards and tested scripts that the
agent would otherwise get wrong, skip, or reinvent badly. Do not restate textbook background
or a library's documentation. If a paragraph would not change the agent's behaviour, cut it.

Write to a smart colleague: explain *why* each rule exists so the agent can generalise to
cases you did not foresee. Avoid shouting (ALWAYS / NEVER in capitals) and rigid scripts
where judgment is needed. Use the imperative mood. No marketing language.

## Folder layout (one skill)

```
skills/<name>/
  SKILL.md            # target 150 to 320 lines, hard cap 500
  references/         # loaded on demand; one topic per file; a '## Contents' list near the top if > 100 lines
    sources.md        # REQUIRED: the bibliography for every claim in this skill (see Evidence)
    <topic>.md
  scripts/            # Python standard library only
    <tool>.py
  evals/              # three test cases in the `claude plugin eval` layout (see Evals)
tests/<name>/         # at the REPO ROOT, not inside the skill: unittest files for the scripts
  test_<tool>.py
```

Reference files are one level deep (SKILL.md links to them; they do not chain to further
files). Every file under `references/` and `scripts/` must be mentioned in SKILL.md with a
sentence saying when to read or run it. No other files (no README inside a skill).

## SKILL.md

Frontmatter (YAML, exactly these keys, in this order):

```yaml
---
name: <folder name: lowercase letters, digits, hyphens>
description: <third person; what it does AND when to use it; concrete trigger terms a scientist would type, including a few Chinese terms; at most 900 characters; no angle brackets>
license: MIT
metadata:
  author: Alex Wang
  version: "0.1.0"
  domain: tumor-immunology
  summary: <one plain sentence for the catalog table, at most 140 characters>
  summary-zh: <the same sentence in Chinese, at most 70 characters>
---
```

Write `description` (and any value containing a colon followed by a space, or ` #`) either as one
double-quoted line or as a `>-` folded block: full YAML parsers reject an unquoted value that
contains `: `, and the skill then fails to load.

Body, in this order (rename headings only if the content truly needs it):

1. `# <Title>` and one short paragraph: the situation, and what the agent must deliver.
2. `## Scope`: what this skill covers, what it does not, and which sibling skill to use
   instead (by name). Three to six lines.
3. `## Workflow`: numbered steps with explicit decision points. Each step says what to do,
   what to check, and why it matters. This is the core of the skill.
4. `## Defaults and decision rules`: compact tables of thresholds and choices, each row with
   a source key such as `[Ayers2017]` that resolves in `references/sources.md`, or the label
   `convention` / `heuristic` when it is not a published threshold.
5. `## Pitfalls`: the mistakes that actually occur (agent mistakes and human mistakes), each
   with the consequence and the fix.
6. `## Scripts`: exact command lines, input format (column names), output format, and how to
   read the result. Scripts are run, not read; say so.
7. `## Report`: what the agent's answer to the researcher must contain (assumptions, data
   checks done, results with uncertainty, limits, what to verify experimentally or
   orthogonally). Give a short template.
8. `## References`: one line per bundled file saying when to open it.

## Evidence rules (the reason this repository exists)

- Every numeric threshold, gene list, antibody clone, dose, reagent choice, guideline
  statement or attributed claim carries a source key. Key format: first author's surname
  (or the issuing body) starting with a capital letter, letters only, then the year, with an
  optional suffix letter: `[Ayers2017]`, `[VanAllen2015]`, `[Jiang2018a]`, `[FDA2020]`.
  In `references/sources.md` each entry is one list item that starts with the bold key, for
  example: `- **[Ayers2017]** Ayers M, et al. IFN-gamma-related mRNA profile predicts clinical
  response to PD-1 blockade. J Clin Invest. 2017;127(8):2930-2940. doi:10.1172/JCI91190. PMID: 28650338.`
  followed, when useful, by one indented line saying what the skill takes from that source. `references/sources.md` lists each key
  with authors, year, title, journal, DOI and (when it exists) PMID.
- You must verify every source you cite against Crossref or PubMed before citing it:
  `curl -s "https://api.crossref.org/works/<DOI>"` (check title, journal, year, first author)
  or `curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=<PMID>&retmode=json"`.
  To find a DOI: `curl -s "https://api.crossref.org/works?query.bibliographic=<title words>&rows=3"`.
  Never cite from memory alone. If you cannot verify a source, drop the claim or mark it
  plainly as `convention (no single primary source)`.
- Verify the *content* of the claim too, not only that the paper exists: open the abstract
  (PubMed efetch, Europe PMC `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=DOI:<doi>&resultType=core&format=json`)
  or the full text when the number matters. A real paper cited for something it does not say
  is the worst kind of error for this repository.
- Separate three kinds of statement and label them: published and validated; common
  convention; this skill's heuristic.
- Date-sensitive facts (regulatory approvals, guideline editions, software versions, database
  releases) carry an "as of <month year>" stamp and an instruction to re-check.
- Gene symbols are current HGNC symbols; give the legacy alias in brackets when people still
  search for it (for example `HAVCR2 [TIM-3]`). State species when it matters (mouse symbols
  in sentence case).
- Research use only. Nothing may read as clinical decision support for an individual patient.
- British or American spelling: use "tumor" in file names, skill names and descriptions (it
  is what people type); either spelling in prose, but be consistent within a file.

## Script rules

- Python >= 3.9, standard library only (no numpy, pandas, scipy). Start with
  `#!/usr/bin/env python3` and a module docstring that includes two usage examples.
- `argparse` with a helpful `--help`; inputs are TSV / CSV / JSON files with documented column
  names; output is JSON or TSV on stdout, or to `--out`.
- Deterministic: no network, no clock, no hidden randomness (any resampling takes `--seed`).
  The single exception is a skill whose purpose is looking things up online: it must declare
  `compatibility: Requires network access to <hosts>` in the frontmatter (between `license` and
  `metadata`), keep all parsing in pure functions, and test those functions on saved responses.
- Fail loudly and usefully: say which column or value is wrong and exit with code 2.
- Put the computation in importable pure functions; keep `main()` thin.
- State in the docstring exactly which formula is implemented and cite its source key.
- Tests (`tests/<name>/test_<tool>.py`, `unittest`): load the script by path with
  `importlib.util`; include at least one golden value derived by hand or taken from a
  published worked example (show the derivation in a comment), edge cases (empty input, ties,
  zeros, a single group), and one command-line smoke test through `subprocess` with
  `sys.executable`. Run them: `python3 -m unittest discover -s tests/<name> -v`.

## Evals (three per skill)

Each case is a folder `skills/<name>/evals/<case-name>/` (lowercase, hyphens, starts with a
letter; never `results`, `mocks`, `files`, `node_modules`, `true`, `false`, `null`):

`prompt.md` (settings block copied exactly, then the prompt as a real user would type it):

```markdown
---
tags: [tumor-immunology]
max_turns: 20
timeout_seconds: 600
allowed_tools: [Read, Glob, Grep, Skill, Write, Edit, Bash]
---

<the user's request, realistic, with concrete details; one of the three cases is in Chinese>
```

`case.yaml`:

```yaml
schema_version: "1.1"
name: <case-name, same as the folder>
expected_outcome: |
  <what a good result looks like, in plain sentences>
```

`graders/NN-few-words.md`, three to six per case, one verifiable statement each:

```markdown
---
type: llm
focus: trace
weight: 1
---

<one objective statement about the result>
```

Graders must discriminate: a capable agent WITHOUT the skill should plausibly fail several
of them (because it lacks the rule, the threshold, or the script), and an agent WITH the
skill should pass. Do not write keyword bingo or statements that any reasonable answer
satisfies. Word each about the result itself, not about where a file is saved.

If a case needs input files, put them in `files/`, end the prompt with a line
`Input files: a.tsv, b.tsv`, add `place-files.sh` beside `prompt.md`:

```bash
#!/usr/bin/env bash
set -eu
cp -R "$(dirname "$0")/files/." .
```

and add to `case.yaml`:

```yaml
context:
  scaffold_script: place-files.sh
```

Keep input files small (a few KB), synthetic, and labelled as synthetic.

## Tone and honesty

- Say what is unknown or contested. Offer the standard options and how to choose.
- Prefer "check X before concluding Y" over rules that pretend biology is tidy.
- Synthetic examples are labelled synthetic. Never present simulated output as a finding.
- No em-dash-heavy prose, no emoji, no hype. Plain, precise scientific English.

## Definition of done (run these from the repository root before you report back)

```bash
python3 tools/validate_skills.py --skill <name>     # zero errors; read the warnings and fix what is real
python3 tools/run_tests.py <name>                   # every test passes
```

Then re-read your SKILL.md once more as the agent that will receive it with no other context:
is every instruction actionable, is anything there only to sound thorough, would you know
what to do at each decision point?
