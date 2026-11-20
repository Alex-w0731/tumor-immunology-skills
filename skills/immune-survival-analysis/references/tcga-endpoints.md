# TCGA endpoint recommendations by cancer type

Condensed from Table 3 of [Liu2018], the TCGA Pan-Cancer Clinical Data Resource (TCGA-CDR).
Each cell gives the authors' recommendation for that endpoint and, in brackets, the number of
events in the TCGA-CDR for that cancer type. "caution" is a recommendation that carries an
asterisk in the paper (read the note); "approx." marks disease-specific survival (DSS) that
had to be approximated because the cause of death was not verified. The resource was
published in April 2018; GDC clinical data may have been updated since, so recount events in
the data you actually analyse (as of November 2026; re-check). Every row below was checked
against Table 3 in the PMC version of the paper.

Endpoints: OS overall survival, PFI progression-free interval, DFI disease-free interval,
DSS disease-specific survival. DFI is defined only for patients documented as disease-free;
de novo stage IV patients and those with less than 90 days of follow-up are excluded from it
[Liu2018].

How to use the table

- Table 3 recommends PFI without reservation for 29 of the 33 types; the exceptions are LAML (no
  PFI), DLBC and KICH (caution) and PCPG (not recommended). The paper's Discussion says 27 of 33,
  which does not match the table. PFI and DFI generally record more events than OS and DSS
  because TCGA follow-up is short [Liu2018]. Prefer PFI when OS is "caution" or
  "no", or when OS events in your analysed subset are too few for the model you plan.
- All four endpoints are recommended without reservation for 13 types: BLCA, CESC, COAD, ESCA,
  HNSC, KIRP, LUAD, LUSC, OV, PAAD, SARC, STAD and UCEC [Liu2018].
- The event counts are for the whole TCGA-CDR cohort. After your exclusions (sample type,
  missing feature, missing covariates) there will be fewer; recompute events per variable.
- The paper's text and Table 3 disagree in a few places: the text does not recommend OS or DSS
  for PRAD while the table marks PRAD OS "caution" (10 events); the text lists KIRC, DLBC and
  PCPG as PFI "with caution" and also flags CHOL and UCS, while the table marks DLBC and KICH
  with an asterisk and PCPG as not recommended. Where they disagree, take the stricter reading
  and say so in the report (heuristic).

| Type | N | OS | PFI | DFI | DSS | Note |
|---|---|---|---|---|---|---|
| ACC | 92 | yes (34) | yes (49) | caution (14) | yes, approx. (30) | few events |
| BLCA | 412 | yes (181) | yes (177) | yes (32) | yes, approx. (124) |  |
| BRCA | 1097 | caution (151) | yes (145) | yes (84) | caution, approx. (83) | OS and DSS need longer follow-up |
| CESC | 307 | yes (71) | yes (71) | yes (26) | yes (54) |  |
| CHOL | 45 | yes (22) | yes (23) | caution (10) | caution, approx. (18) | small sample for every endpoint |
| COAD | 459 | yes (102) | yes (123) | yes (24) | yes, approx. (64) |  |
| DLBC | 48 | no (9) | caution (12) | no (4) | no (4) | small sample, few events, short follow-up |
| ESCA | 185 | yes (77) | yes (87) | yes (23) | yes, approx. (51) |  |
| GBM | 596 | yes (491) | yes (506) | no (2) | yes, approx. (445) | few disease-free cases (DFI) |
| HNSC | 528 | yes (223) | yes (198) | yes (28) | yes, approx. (130) |  |
| KICH | 113 | caution (13) | caution (17) | no (6) | caution, approx. (10) | very few events, short follow-up |
| KIRC | 537 | yes (177) | yes (162) | caution (15) | yes, approx. (110) | few events |
| KIRP | 291 | yes (44) | yes (58) | yes (28) | yes, approx. (28) |  |
| LAML | 200 | yes (133) | n/a | n/a | n/a | OS only |
| LGG | 515 | caution (125) | yes (192) | yes (20) | caution, approx. (113) | OS and DSS need longer follow-up |
| LIHC | 377 | yes (132) | yes (185) | yes (147) | caution, approx. (80) | DSS needs longer follow-up |
| LUAD | 522 | yes (188) | yes (213) | yes (92) | yes, approx. (116) |  |
| LUSC | 504 | yes (219) | yes (149) | yes (63) | yes, approx. (91) |  |
| MESO | 87 | yes (74) | yes (61) | no (7) | yes, approx. (43) | small sample for DFI |
| OV | 587 | yes (349) | yes (414) | yes (196) | yes, approx. (302) |  |
| PAAD | 185 | yes (100) | yes (110) | yes (23) | yes (79) |  |
| PCPG | 179 | no (6) | no (21) | no (4) | no (4) | every endpoint needs longer follow-up; few events |
| PRAD | 500 | caution (10) | yes (93) | yes (30) | no (5) | OS and DSS need longer follow-up |
| READ | 170 | caution (26) | yes (39) | no (7) | caution, approx. (15) | OS, DSS and DFI need longer follow-up; very few DFI events |
| SARC | 261 | yes (99) | yes (139) | yes (67) | yes, approx. (81) |  |
| SKCM | 470 | yes (216) | yes (309) | n/a | yes, approx. (190) | DFI cannot be derived |
| STAD | 443 | yes (172) | yes (143) | yes (46) | yes, approx. (103) |  |
| TGCT | 134 | no (4) | yes (35) | yes (27) | no (3) | few OS and DSS events; short follow-up |
| THCA | 507 | caution (16) | yes (52) | yes (26) | no (7) | few OS and DSS events; short follow-up |
| THYM | 124 | no (9) | yes (22) | n/a | no (4) | very few OS and DSS events; short follow-up; DFI cannot be derived |
| UCEC | 548 | yes (91) | yes (124) | yes (57) | yes, approx. (60) |  |
| UCS | 57 | yes (35) | yes (37) | caution (10) | yes, approx. (31) | small sample |
| UVM | 80 | yes (23) | yes (30) | n/a | yes (21) | DFI cannot be derived |
