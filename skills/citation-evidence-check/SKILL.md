---
name: citation-evidence-check
description: >-
  Verifies the references and cited claims in tumor immunology manuscripts, grant text, slides
  and agent-written reports before they are shared: each reference resolves at Crossref or
  PubMed, its title, year, first author and journal match, retraction, expression-of-concern and
  correction notices are found, and the statement attributed to it is supported by what the paper
  reports, with numbers, species and model system, and primary data versus review or preprint
  checked. Ships a tested lookup script with saved-response replay. Use when the user asks to
  check references, verify citations, find fabricated or hallucinated references, check for
  retracted papers, fact-check a draft or an AI-written text, or confirm that a paper says what
  the text claims; also for 核对参考文献, 引用核查, 文献真实性, 撤稿, AI编造文献, 事实核查.
license: MIT
compatibility: Requires network access to api.crossref.org, doi.org, eutils.ncbi.nlm.nih.gov and www.ebi.ac.uk
metadata:
  author: Alex Wang
  version: "0.1.1"
  domain: tumor-immunology
  summary: Checks that references exist, match their metadata, are not retracted, and actually support the claims they are cited for.
  summary-zh: 核对参考文献是否真实、信息无误、未撤稿，并确认原文确实支持所引用的论述。
---

# Citation and evidence check

A tumour-immunology text that cites the literature (manuscript, grant section, slides, or a
report written by an agent) is about to be shared. Deliver, for every citation, two separate
answers: is the reference real and correctly described (it resolves, its details match, it is
not retracted), and does the paper support the statement attached to it (checked against the
abstract, or the full text when a number or detail depends on it). Then propose the smallest
corrections and list what could not be verified. ICMJE expects authors to be able to attest that
each reference supports the statement it is attached to [ICMJE2026]; this check is how they can.
Never repair or confirm a citation from memory: language models, including you, produce
plausible references that do not exist [Walters2023].

## Scope

- Covers: reference lists in any citation style, prose, slide text and agent output; existence
  and metadata, retraction, expression-of-concern and correction notices; evidence type
  (primary, review, preprint, editorial, trial); claim support including numbers, population,
  species and model system; statements that carry no citation.
- Not covered: appraising a cited analysis in depth (use `immune-survival-analysis` for hazard
  ratios and survival claims, `response-biomarker-evaluation` for response-prediction claims,
  `single-cell-patient-level-stats` for single-cell comparisons, `immune-signature-curation`
  for gene-signature provenance); bibliography formatting; plagiarism checks.

## Workflow

1. **Inventory statements and references before judging anything.** Make one row per
   (statement, citation) pair with `id`, the statement verbatim in `claim`, and whatever the text
   gives: `doi`, `pmid`, `title`, `year`, `first_author`, `journal`
   (`references/claim-review.md` section 1 maps BibTeX, RIS and Vancouver fields). Keep the
   statement exact; a paraphrase hides the overstatement you are looking for. Add rows without a
   reference for uncited statements that carry a number, a causal verb ("drives", "determines",
   "is required for") or a general claim about patients. Why: citations fail in two independent
   ways, as wrong references and as real references that do not say what is attributed to them;
   in medical journals, estimates of quotation errors range from 14.5% of quotations [Mogull2017]
   to 25.4%, 11.9% major [Jergas2015], depending on how errors are counted.
   - If a language model wrote the text, treat every reference as unverified however plausible
     it looks. In one study 55% of GPT-3.5 and 18% of GPT-4 citations were fabricated
     [Walters2023]; in another, 47% of references were fabricated, 46% real but inaccurate, and a
     wrong PMID was the commonest error [Bhattacharyya2023]; the authors of a third evaluation
     conclude that model-generated references need thorough validation [Chelli2024].

2. **Run the lookup** (command in Scripts) with `--abstracts --save responses/`. The saved folder
   is the provenance record; `--offline responses/` reproduces the report without network.
   Only reference metadata (title, first author, year, DOI, PMID) leaves the machine; the `claim`
   column is never sent. Ask first if the reference list itself is confidential.
   - Decision: references with `lookup_failed` are unchecked. Rerun later; never fill the gap
     from memory or a search-engine snippet, and never report them as verified.
   - Decision: `unconfirmed` means nothing was compared (identifier only, or a DOI registered
     outside Crossref). Compare the registered title and first author with the text yourself
     (`record`, or the doi.org landing page) before treating the reference as identified.

3. **Resolve identity problems** with the verdict table below. Correct a reference only to the
   work the text demonstrably means: same title, first author and year. When the cited work
   cannot be found after the searches in `references/claim-review.md` section 2, report it as
   "not found: probably fabricated or garbled". Do not substitute the closest real paper; the
   script refuses loose title matches for the same reason. If another paper supports the
   statement, offer it as a new citation, labelled new, after checking its content (step 5).

4. **Act on notices** (`notice_status`, table below). Retracted papers keep being cited: in
   biomedicine only 5.4% of post-retraction citation contexts acknowledged the retraction
   [Hsiao2021]. Authors are responsible for not citing retracted work except when discussing the
   retraction, and ICMJE treats PubMed as the authoritative source on retractions in
   MEDLINE-indexed journals [ICMJE2026]; the script reads PubMed's publication types and links
   [NLM2026] and Crossref's publisher and Retraction Watch entries [Crossref2026b].
   `none_reported` means neither database reported a notice on the day of the check; it says
   nothing about whether the paper is sound. `unknown` means the cited work was not identified,
   so settle identity (step 3) first.

5. **Check that each paper supports its statement.** Read the abstract the script placed next to
   the claim. Open the full text (PMC when `pmcid` is present, otherwise the publisher) when the
   statement depends on a number, a subgroup, a figure, a method, or anything the abstract does
   not state. Compare, in this order:
   - numbers: value, denominator, and which cohort, arm, tumour type, dose or time point they
     belong to. `claim_numbers_not_in_abstract` lists claim numbers absent from the abstract; a
     number that is present may still belong to another subgroup;
   - direction and size of the effect;
   - population: cancer type, stage, line of therapy, treatment;
   - species and model system (step 7);
   - design: an association does not support a causal verb;
   - ownership: if the paper only cites the finding (introduction, discussion), it is a secondary
     citation and the statement needs the primary source.
   Read the results, not only the concluding sentence: reports can present findings more
   favourably than the data allow ("spin") [Boutron2018]. Classify each pair with the support
   table below. Saying that a paper does not report something (a rate, a subgroup, a result)
   needs the full text; from the abstract alone write "not in the abstract", because an abstract
   carries only part of the results.

6. **Check the evidence type** from `flags` and the paper itself.
   - Review: acceptable for a concept or an overview. For a specific finding, trace it to the
     primary paper the review cites and check it there [ICMJE2026]; chains of citation can turn a
     hypothesis into apparent fact [Greenberg2009].
   - Preprint: label it as a preprint [ICMJE2026]. With `published_version_available`, check the
     statement against the published version and cite that; abstract conclusions changed in 7.2%
     of non-COVID-19 preprints that were later published, mostly not qualitatively [Brierley2022].
   - Editorial, comment or letter: opinion, not evidence of an effect. Case report: one instance,
     not a rate. Clinical trial: check the phase, arm and endpoint behind the cited number. Trial
     protocol (`trial_protocol`): a design, not results, so it cannot support an outcome.

7. **State species and model system for every generalised statement.** Tumour-immunology claims
   are often built on mouse data and written about patients. Record what the supporting paper
   used: patients (cohort, cancer type, n), human tissue or cells ex vivo, syngeneic mouse tumour
   lines (MC38, B16, CT26, 4T1), genetically engineered mice, xenografts or PDX, or in vitro
   systems. `model_cues` lists species and model words in the abstract; confirm in the methods.
   A xenograft in an immunodeficient host cannot show a T-cell-mediated effect unless human
   immune cells were supplied (humanised mice, transferred T cells). Mouse and human immunology
   differ in leukocyte subsets, receptors, cytokines and costimulatory molecules, which matters
   whenever mice model human disease [Mestas2004]. In one comparison, language-model summaries of
   research had nearly five times the odds of containing broad generalisations that human-written
   summaries had (OR 4.85) [Peters2025], so check agent-written text sentence by sentence. Fix with
   "in mouse models, ..." or a human source.

8. **Decide what happens to anything unverified**, and say which option you chose: narrow the
   statement to what a verified source shows; mark it for the author (`[unverified: reason]`);
   delete it; or propose a replacement source you looked up and read. Never keep an unverified
   statement silently and never cite from memory. When only the abstract was accessible, say so
   for that statement.
   - Keep the corrected text to the author's own statements, fixed, narrowed, marked or removed.
     List replacement sources and any new sentence beside it, marked new, with the identifier and
     what you checked it against (abstract or full text), unless the user asked you to insert
     replacements. A new paper inside a corrected paragraph reads as the author's checked claim,
     and co-authors cannot tell it from the original.

9. **Report** with the template below. Show proposed changes against the original wording; in
   someone else's text the author decides.

## Defaults and decision rules

| Verdict | Meaning | Action |
|---|---|---|
| verified | identifiers resolve; title agreement at least 0.9 (when a title was given); year within one; author and journal agree or were not given | go on to claim support |
| check_title | title agreement 0.6 to 0.9 | compare the titles: a subtitle, a translated title, or a different paper from the same group |
| metadata_errors | year differs by more than one, or first author or journal disagree | correct from the registered record after confirming it is the intended work |
| wrong_identifier | the cited work exists, but a DOI or PMID is broken or belongs to another work | replace the identifier; check the claim against the real work |
| identifier_conflict | DOI and PMID point to different works | decide from title and authors which work is meant; check that one |
| mismatch | the identifier resolves to an unrelated work and the Crossref title search found no matching record; `notice_status` is then `unknown` | treat as fabricated or garbled; search (claim-review.md section 2); never substitute |
| unconfirmed | only an identifier was given, or the DOI is registered outside Crossref (DataCite, for example) and nothing was compared | compare the registered title and authors with what the text describes |
| not_found | nothing resolves and the title search found nothing close | as mismatch |
| lookup_failed | service error or missing saved response | unchecked; rerun |

| Rule | Value | Source |
|---|---|---|
| Title agreement (Dice on title words) | match at least 0.9; partial 0.6 to 0.9; different work below 0.6 | heuristic; coefficient [Dice1945] |
| Title-search hit accepted as the cited work | agreement at least 0.9, or 0.6 to 0.9 with year within one and the same first author | heuristic |
| Year off by one | not an error (online and print years differ) | convention |
| Retracted | no support for any claim; cite only when discussing the retraction, with the notice | [ICMJE2026] |
| Expression of concern | not sole support; tell the author | convention |
| Correction | read it; it matters only if it changes the cited result | convention |
| Cited item is itself a notice (`is_notice`) | cite the article it concerns for findings; cite the notice only when writing about the retraction or correction | convention |
| Review cited for a specific finding | trace to the primary source | [ICMJE2026] |
| Preprint | label it; cite the published version once the statement is checked there | [ICMJE2026] [Brierley2022] |
| Request pacing | 0.4 s between requests (NCBI: three per second without an API key, as announced in 2017; re-check the E-utilities documentation); back off on HTTP 429 | [NCBI2017] [Crossref2026a] |

| Support category | Definition | Action |
|---|---|---|
| supported | the paper's own results state it, for the same population, species and design | keep |
| overstated | oversimplified or overgeneralised, Mogull's "minor" error [Mogull2017]. This skill files here: a wider species, population or cancer type, a causal verb for an association, and a rate generalised from a few cases | narrow the wording; still a required fix |
| not supported | the source fails to substantiate the statement, is unrelated to it, or contradicts it: Mogull's "major" error [Mogull2017]; also a number the paper does not report | remove, or replace with a checked source |
| secondary citation | the paper states it only by citing other work (introduction, discussion, review); Mogull counts improper secondary citation separately [Mogull2017] | find the primary source, check the statement there, cite that |
| not assessable | abstract silent and full text not accessible | report "abstract only, not assessable"; never count as supported |

## Pitfalls

- **Existence taken as support.** A resolving DOI says nothing about the statement, and a real
  paper cited for something it does not say is the error readers are least likely to catch.
  Keep identity and support as separate columns.
- **Repairing an invented reference.** Matching it to the nearest real title, author or DOI
  attaches a real paper to a statement nobody checked. Keep "not found" unless title, first
  author and year agree; anything else is a new citation that must be read first.
- **Verifying from memory.** Recall is where fabricated references come from. Every verdict rests
  on a lookup result or a document read in this session.
- **Number found, wrong context.** The digits are in the abstract but for another tumour type,
  arm, dose or time point. Compare the sentence around the number.
- **Small denominators turned into rates.** A percentage built from a few patients is not a
  population estimate; give the counts.
- **Mouse result written about patients**, or an immune mechanism claimed from xenografts in
  immunodeficient mice. Apply step 7.
- **Retraction checked in one database.** Crossref and PubMed are maintained separately and need
  not agree on a given day; the script queries both, and `none_reported` covers those two only.
- **Citing a review, a preprint or a citing paper for a primary result.** Trace the claim to the
  paper that produced the data.
- **Over-flagging formatting.** Greek letters (IFN-γ versus IFN-gamma), journal abbreviations
  and one-year differences are not errors; the script normalises them. Noise buries real problems.
- **Silent rewriting.** Changing an author's sentence without showing the original and the reason
  removes their chance to disagree.

## Scripts

Run the script (path relative to this skill's folder); there is no need to read its source.

```bash
python3 scripts/verify_references.py refs.csv --abstracts --save responses/ --out report.json
python3 scripts/verify_references.py refs.csv --offline responses/ --format tsv --out report.tsv
```

- Input: CSV, TSV or JSON (a list of objects, or `{"references": [...]}`), UTF-8. Columns, case
  insensitive: `doi`, `pmid`, `title`, `year`, `first_author` (or `author`, or `authors` from
  which the first name is taken), optional `id`, `journal`, `claim`. DOIs may be URLs or carry a
  `doi:` prefix; a year may carry a suffix or a date (`2017a`, `2017 Jun 26`). A malformed DOI,
  PMID or year is reported in `issues` and ignored. The run stops with exit code 2, before any
  lookup, when a row has none of DOI, PMID and title, when the file has no doi, pmid or title
  column or is not UTF-8, or when an option or the `--out` folder is invalid.
- Output (JSON): `summary` with counts per verdict and notice status, then per reference
  `verdict`; `notice_status` (retracted, expression_of_concern, is_notice, corrected,
  none_reported, or unknown when the cited work was not identified); `notices` (kind, source,
  notice identifier, date); `notice_for` (when the cited item is itself a notice); `flags`
  (preprint, published_version_available, review, editorial_or_comment, clinical_trial,
  trial_protocol, case_report, title_says_retracted, title_says_notice,
  not_a_journal_article:TYPE); `title_score`, `year_check`, `author_check`, `journal_check`;
  `record` (registered DOI, PMID, title, years, first author, types); `issues` in plain
  language; with `--abstracts` also `abstract`, `pmcid`, `open_access`, `model_cues` and
  `claim_numbers_not_in_abstract`. `--format tsv` gives one row per reference without abstracts.
- Options: `--delay` seconds between requests (default 0.4), `--timeout`. Exit code 0 whenever
  the run completes, whatever the verdicts, so read `summary` rather than the exit code.
- Reading it: the script decides identity and notices, never support. `model_cues` and
  `claim_numbers_not_in_abstract` are prompts for your reading: numbers written as words are not
  compared, both ends of a range (26-73%) are, and the abstract may be missing (`abstract`
  empty), in which case open the paper. The mouse cues are a hand-made list of common lines and
  strains; a missing cue does not mean human data.
- The Crossref and PubMed notice fields were checked in November 2026. If a paper you know to be
  retracted comes back `none_reported`, check its PubMed and publisher pages by hand and re-test
  the script on that DOI: the services may have changed their fields.

## Report

State which version of the text was checked and which references you converted by hand. Lead
with what must change, because that is what the researcher acts on. Then:

```markdown
## Reference and claim check: <document>, <date>
Checked N references and M statements (K without a citation) against Crossref, PubMed and
Europe PMC on <date>; raw responses saved in <folder>.

Must change before sharing: [n] retracted; [n] not found, probably fabricated; [n] not supported.
Should change: [n] overstated (species, population, causal wording); [n] metadata errors.

| # | Statement (short) | Reference as cited | Identity | Notices | Support | Evidence type; species/model | Checked against | Action |
|---|---|---|---|---|---|---|---|---|
| 1 | ... | Author Year, journal | verified | none reported | overstated | primary; mouse (MC38) | abstract | narrow to mouse models |

### Proposed changes
- [n] original -> proposed. Reason.

### Optional new sources (not in the corrected text)
- [n] new: Author Year, DOI or PMID; checked against abstract or full text; what it supports.

### Unresolved, for the author
- Statements that could not be verified, what was searched, and why it failed.

### Limits
Numbers checked against abstracts only where stated; sources without identifiers (books,
guidelines) checked by hand or not at all; "none reported" means none in Crossref or PubMed on
<date>.
```

## References

- `references/claim-review.md`: open when converting a reference list to script input, when a
  reference was not found and needs a second search, when the text cites GEO accessions or trial
  numbers, and when judging support (failure patterns and wording fixes).
- `references/sources.md`: the bibliography behind every rule here, with what was verified and
  how far each source was read.
