# Claim review: inputs, searches and reading for support

## Contents

1. Turning a reference list into script input
2. Searching for a reference the script did not find
3. Identifiers that are not articles (datasets, trials)
4. Reading a paper for support: what to compare
5. Wording fixes that keep the author's point

## 1. Turning a reference list into script input

Copy fields as they are written; do not correct anything while converting, or the check will
confirm your correction instead of the author's text.

| Source format | doi | pmid | title | year | first_author | journal |
|---|---|---|---|---|---|---|
| BibTeX | `doi` | `pmid` or a PMID in `note` | `title` (drop braces) | `year` | first name in `author` | `journal` |
| RIS | `DO` | `AN` only if it is a PMID | `TI` or `T1` | `PY` or `Y1` (first 4 digits) | first `AU` | `JO`, `JF` or `T2` |
| Vancouver or APA string | after "doi:" or doi.org/ | after "PMID:" | the sentence after the authors | the 4-digit year | the first surname with initials | the abbreviation before the year |

Put the sentence that cites the reference in `claim`, verbatim. A reference cited by several
sentences gets one row per sentence, with the same `id` plus a suffix (`3a`, `3b`).

## 2. Searching for a reference the script did not find

The script already ran a Crossref bibliographic search on title, first author and year. Before
calling a reference fabricated, try two independent searches (syntax confirmed November 2026).
With curl, pass `-g` (or write the brackets as `%5B` and `%5D`): otherwise curl reads `[1au]` as
a URL pattern and returns nothing, which looks like "no match".

- PubMed (E-utilities [NCBI2026]): `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&term=Surname[1au]+AND+2016[dp]+AND+keyword[ti]`
  (`[1au]` first author, `[dp]` publication date, `[ti]` title word; widen by dropping one term).
- Europe PMC (also indexes preprints): `https://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&resultType=lite&query=AUTH:"Surname AB" AND PUB_YEAR:2016 AND TITLE:"keyword"`;
  add `AND SRC:PPR` to search preprints only.
- Abstract as plain text: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=PMID&rettype=abstract&retmode=text`.

A hit is the cited work only if title, first author and year all agree (one-year difference
allowed). A paper on the same topic by other authors is not the cited work, even when it would
support the sentence; offer it as a new citation after reading it. Record what you searched, so
the author can repeat it: "PubMed and Europe PMC searched for first author, year and title words,
no match" is a finding; "could not find" is not.

## 3. Identifiers that are not articles (datasets, trials)

Tumour-immunology text often cites accessions and trial numbers. The script does not check these;
look them up directly (endpoints confirmed November 2026):

- GEO series: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmode=json&term=GSE72056[ACCN]+AND+gse[ETYP]`,
  then `esummary.fcgi?db=gds&retmode=json&id=UID` for title, organism (`taxon`), sample count and
  linked PMIDs. Check that organism and sample count match what the text says.
- ClinicalTrials.gov: `https://clinicaltrials.gov/api/v2/studies/NCT00730639?fields=NCTId,BriefTitle,OverallStatus,Phase`.
  Check that the phase, condition and intervention match the sentence that cites the trial.

## 4. Reading a paper for support: what to compare

Work from the paper's own results. Typical failures, with what to look for:

| Failure | What to look for in the source | Category |
|---|---|---|
| Number from another subgroup | The sentence around the number: tumour type, arm, dose, time point, denominator | not supported if the number belongs elsewhere |
| Rate from a handful of cases | Sample size behind a percentage; "three of eight patients" is not a population rate | overstated if the counts are the paper's (give them); not supported if the paper reports no such rate |
| Association written as cause | Design: observational cohort, correlation, signature association | overstated |
| Mouse result written about patients | Species in methods; syngeneic lines (MC38, B16, CT26, 4T1), GEMMs, germ-free or antibiotic-treated mice | overstated |
| Immune mechanism from xenografts | Immunodeficient host without human immune reconstitution cannot show a T-cell-mediated effect | overstated or not supported |
| In vitro result written as in vivo | Cell lines, co-culture, organoids without immune cells | overstated |
| Finding the paper only cites | Statement appears in the introduction or discussion with its own citation | secondary citation (check and cite the primary source) |
| Conclusion stronger than results | Abstract conclusion versus the reported effect sizes and confidence intervals | overstated |
| Single cancer type generalised | "in melanoma" in the source, "in solid tumours" in the text | overstated |

When the abstract is silent and the full text is not accessible, the category is "not
assessable", not "supported".

## 5. Wording fixes that keep the author's point

Narrow rather than delete when the evidence supports a smaller claim:

- Add the model: "In syngeneic mouse models, treatment with X improved tumour control [n]."
- Replace a causal verb: "is associated with" for observational data; keep "drives" only for
  perturbation data in the stated system.
- Give the denominator: "in 3 of 12 evaluable patients [n]".
- Restrict the population: "in non-small-cell lung cancer" instead of "in solid tumours".
- Mark what remains open for the author: `[unverified: no source found for this rate]`.

Each fix goes in the report as original, proposed, reason. The author decides.
