# Response labels for biomarker evaluation

Open this file when turning a dataset's response column into responder / non-responder
labels, or when two cohorts label response differently.

## Mapping by source field

| Source field | Responder (1) | Non-responder (0) | Set aside, count, test in a sensitivity analysis | Source |
|---|---|---|---|---|
| RECIST 1.1 best overall response | CR, PR | SD, PD | NE, not assessed, no post-baseline scan | [Eisenhauer2009] |
| iRECIST best overall response (iBOR) | iCR, iPR | iSD, iUPD, iCPD | NE | [Seymour2017] |
| Clinical-benefit label | as the source study defines it, e.g. PR or SD lasting more than 6 months | e.g. no durable benefit | follow-up shorter than the benefit window without progression | [Rizvi2015] |
| Binary column already supplied | read the paper's definition first | | missing values | convention |

Notes on the table:

- Confirmation. RECIST 1.1 requires responses to be confirmed when response is the primary
  endpoint of a non-randomized trial, but not in randomized trials [Eisenhauer2009]. A pooled
  dataset can therefore mix confirmed and unconfirmed responses. If unconfirmed responses are
  flagged (for example uPR), keep them as their own category and run the analysis with them
  counted both ways (convention).
- iRECIST and pseudo-progression. Under iRECIST an iUPD does not override a later iSD, iPR or
  iCR in the best overall response [Seymour2017]; under RECIST 1.1 the same patient is PD. A
  pseudo-progressor can be a responder in an iRECIST cohort and a non-responder in a RECIST
  1.1 cohort, so record which criteria each cohort used before pooling or comparing.
- Clinical-benefit windows differ between studies. Take each study's own definition, write it
  down, and never relabel patients with too little follow-up as non-responders: Rizvi et al.
  kept patients who had not reached 6 months of follow-up apart from both groups [Rizvi2015].
- The tmelens IMvigor210 example codes `binaryResponse == "CR/PR"` as 1, SD/PD as 0 and a
  missing value (NE) as excluded. Column names in public releases change; check them after
  every download, as tmelens advises.

## Sensitivity-analysis plan

Run each row with `scripts/roc_auc.py`, changing only `--positive`, `--negative` or the input
rows, and report every AUC with its interval.

| Analysis | Responders | Non-responders | Excluded | What it shows |
|---|---|---|---|---|
| Primary | CR, PR | SD, PD | NE | the main estimate |
| S1 | CR, PR | SD, PD, NE | none | NE treated as failures (intention-to-treat style) |
| S2 | clinical benefit | no clinical benefit | follow-up too short | whether durable disease control tells the same story |
| S3, descriptive only | CR, PR | PD | SD, NE | contrast of extremes; expected to overstate [Lijmer1999] |
| S4 | as primary | as primary | none by timing | effect of including on-treatment samples, if any exist |

```bash
python3 scripts/roc_auc.py cohort.tsv --score s --label bor --positive CR,PR --negative SD,PD --seed 1      # primary
python3 scripts/roc_auc.py cohort.tsv --score s --label bor --positive CR,PR --negative SD,PD,NE --seed 1   # S1
python3 scripts/roc_auc.py cohort.tsv --score s --label bor --positive CR,PR --negative PD --seed 1         # S3
```

If the AUC moves by more than the half-width of its own interval between the primary analysis
and S1 or S2, the label definition is driving the result; say so (heuristic).

## Sample timing and other checks

- Baseline samples only in the primary analysis. Record the time from biopsy to the first dose
  and any treatment given in between: an archival biopsy taken before an intervening line of
  therapy may not describe the tumour that met the checkpoint inhibitor (convention).
- On-treatment and post-progression samples answer a different question (what changes with
  response) and can reflect the response itself; analyse them separately.
- One patient, one row. If a patient contributes several lesions or time points, choose the
  prespecified one before looking at outcomes.
- Report the counts behind every label category, including the excluded ones. STARD 2015
  asks how indeterminate results and missing data were handled (items 15 and 16) and for a
  flow diagram of participants (item 19) [Bossuyt2015].
