---
type: llm
focus: trace
weight: 1
---

The cleaned GMT keeps the two spreadsheet-damaged APC_MHCII fields (8-Mar and 43891) out of the set until the original source is checked: its APC_MHCII line contains neither the damaged values nor symbols inferred from them (MARCHF8, MARCHF1, MARCH8, MARCH1 or another candidate), and the reply lists both fields as lookups to make before scoring. Writing an inferred symbol into the cleaned GMT fails, even when the reply or a change log asks the user to confirm it.
