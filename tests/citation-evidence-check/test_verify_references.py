"""Unit tests for skills/citation-evidence-check/scripts/verify_references.py.

No network is used. Responses below are trimmed copies of live Crossref, PubMed and Europe PMC
records retrieved in November 2026 (only the fields the script reads are kept), plus a few
records marked SYNTHETIC where a structure had to be shown that no convenient real record has.
"""
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
SCRIPT = os.path.join(ROOT, "skills", "citation-evidence-check", "scripts", "verify_references.py")

_spec = importlib.util.spec_from_file_location("verify_references", SCRIPT)
vr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(vr)

# --- trimmed live records (November 2026) ---------------------------------------------------
CR_AYERS = {"DOI": "10.1172/jci91190", "type": "journal-article",
            "title": ["IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade"],
            "container-title": ["Journal of Clinical Investigation"],
            "author": [{"family": "Ayers", "sequence": "first"}, {"family": "Ribas", "sequence": "additional"}],
            "issued": {"date-parts": [[2017]]}}
CR_KURTULUS = {"DOI": "10.1172/jci81187", "type": "journal-article",
               "title": ["TIGIT predominantly regulates the immune response via regulatory T cells"],
               "container-title": ["Journal of Clinical Investigation"],
               "author": [{"family": "Kurtulus", "sequence": "first"}], "issued": {"date-parts": [[2015]]},
               "updated-by": [{"DOI": "10.1172/jci183278", "type": "retraction", "label": "Retraction",
                               "source": "retraction-watch", "updated": {"date-parts": [[2024, 7, 1]]}, "record-id": "54778"}]}
CR_GROOM = {"DOI": "10.1038/s41590-021-00999-x", "type": "journal-article",
            "title": ["Diversity in science requires mentoring for all, by all"], "container-title": ["Nature Immunology"],
            "author": [{"family": "Groom", "sequence": "first"}], "issued": {"date-parts": [[2021]]}}
CR_IMMUNITY = {"DOI": "10.1016/j.immuni.2013.11.005", "type": "journal-article",
               "title": ["RETRACTED: FoxP3 Acts as a Cotranscription Factor with STAT3 in Tumor-Induced Regulatory T Cells"],
               "author": [{"family": "Hossain", "sequence": "first"}], "issued": {"date-parts": [[2013]]},
               "update-to": [{"DOI": "10.1016/j.immuni.2013.11.005", "type": "retraction", "source": "publisher",
                              "updated": {"date-parts": [[2013, 12, 12]]}}],
               "updated-by": [{"DOI": "10.1016/j.immuni.2021.08.008", "type": "retraction", "source": "retraction-watch",
                               "updated": {"date-parts": [[2021, 9, 14]]}},
                              {"DOI": "10.1016/j.immuni.2013.11.005", "type": "retraction", "source": "publisher",
                               "updated": {"date-parts": [[2013, 12, 12]]}}]}
CR_PREPRINT = {"DOI": "10.2196/preprints.53164", "type": "posted-content", "subtype": "preprint",
               "title": ["Hallucination Rates and Reference Accuracy of ChatGPT and Bard for Systematic Reviews: "
                         "Comparative Analysis (Preprint)"],
               "author": [{"family": "Chelli", "sequence": "first"}], "issued": {"date-parts": [[2023]]},
               "relation": {"is-preprint-of": [{"id-type": "doi", "id": "10.2196/53164", "asserted-by": "subject"}]}}
PM_AYERS = {"uid": "28650338", "title": "IFN-γ-related mRNA profile predicts clinical response to PD-1 blockade.",
            "pubdate": "2017", "authors": [{"name": "Ayers M", "authtype": "Author"}],
            "pubtype": ["Clinical Trial, Phase I", "Journal Article", "Randomized Controlled Trial"],
            "articleids": [{"idtype": "pubmed", "value": "28650338"}, {"idtype": "doi", "value": "10.1172/JCI91190"}]}
PM_KURTULUS = {"uid": "26413872", "title": "TIGIT predominantly regulates the immune response via regulatory T cells.",
               "pubdate": "2015 Nov 2", "authors": [{"name": "Kurtulus S", "authtype": "Author"}],
               "pubtype": ["Journal Article", "Research Support, N.I.H., Extramural", "Retracted Publication"],
               "articleids": [{"idtype": "doi", "value": "10.1172/JCI81187"}],
               "references": [{"refsource": "J Clin Invest. 2024 Jul 1;134(13):e183278. doi: 10.1172/JCI183278.",
                               "reftype": "Retraction in", "pmid": 38949028, "note": ""}]}
PM_CHEN = {"uid": "28102259", "title": "Elements of cancer immunity and the cancer-immune set point.", "pubdate": "2017",
           "authors": [{"name": "Chen DS", "authtype": "Author"}, {"name": "Mellman I", "authtype": "Author"}],
           "pubtype": ["Journal Article", "Review"], "articleids": [{"idtype": "doi", "value": "10.1038/nature21349"}]}
PM_DUHEN = {"uid": "30006565", "title": "Co-expression of CD39 and CD103 identifies tumor-reactive CD8 T cells in human solid tumors.",
            "pubdate": "2018", "source": "Nat Commun", "fulljournalname": "Nature communications",
            "authors": [{"name": "Duhen T", "authtype": "Author"}], "pubtype": ["Journal Article"],
            "articleids": [{"idtype": "doi", "value": "10.1038/s41467-018-05072-0"}]}
# Abstract excerpt of Topalian et al., N Engl J Med 2012 (PMID 22658127), as served by Europe PMC.
TOPALIAN_ABSTRACT = ("<h4>Results</h4>A total of 296 patients received treatment through February 24, 2012. "
                     "Cumulative response rates (all doses) were 18% among patients with non-small-cell lung cancer "
                     "(14 of 76 patients), 28% among patients with melanoma (26 of 94 patients)")
EPMC_TOPALIAN = {"hitCount": 1, "resultList": {"result": [{"pmid": "22658127", "pmcid": "PMC3544539",
                                                          "abstractText": TOPALIAN_ABSTRACT,
                                                          "pubTypeList": {"pubType": ["Journal Article"]}}]}}
# --- SYNTHETIC records (structure only) -------------------------------------------------------
CR_NOTICE = {"DOI": "10.5555/synthetic.notice", "type": "journal-article", "title": ["Retraction notice to: An example article"],
             "update-to": [{"DOI": "10.5555/synthetic.article", "type": "retraction", "source": "publisher",
                            "updated": {"date-parts": [[2020, 1, 2]]}}]}
CR_CORRECTED = {"DOI": "10.5555/synthetic.corrected", "type": "journal-article", "title": ["An example article"],
                "updated-by": [{"DOI": "10.5555/synthetic.erratum", "type": "correction", "source": "publisher",
                                "updated": {"date-parts": [[2021, 3, 4]]}}]}


def ref(**kw):
    base = {"id": "1", "doi": "", "pmid": "", "title": "", "year": None, "first_author": "", "journal": "",
            "claim": "", "input_issues": []}
    base.update(kw)
    return base


def write(folder, name, payload):
    path = os.path.join(folder, name)
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(payload if isinstance(payload, str) else json.dumps(payload))
    return path


def save(folder, kind, key, payload):
    write(folder, vr.cache_name(kind, key), payload)


def esearch(*ids):
    return {"esearchresult": {"count": str(len(ids)), "idlist": list(ids)}}


def esummary(rec):
    return {"result": {"uids": [rec["uid"]], rec["uid"]: rec}}


class TitleAgreement(unittest.TestCase):
    def test_golden_dice_with_greek_letter(self):
        # A = {ifn, gamma, related, mrna, profile, predicts, clinical, response, pd, 1, blockade} (11 words; 'to' dropped)
        # B = A without 'clinical' (10 words); |A & B| = 10, Dice = 2*10/(11+10) = 20/21
        score = vr.dice("IFN-γ–related mRNA profile predicts clinical response to PD-1 blockade",
                        "IFN-gamma related mRNA profile predicts response to PD-1 blockade")
        self.assertAlmostEqual(score, 20 / 21)

    def test_golden_partial_dice(self):
        # A = {pd, 1, blockade, induces, responses, inhibiting, adaptive, immune, resistance} (9; 'by' dropped)
        # B = {pd, 1, blockade, induces, durable, responses, melanoma} (7; 'in' dropped); shared 5 -> 10/16
        score = vr.dice("PD-1 blockade induces responses by inhibiting adaptive immune resistance",
                        "PD-1 blockade induces durable responses in melanoma")
        self.assertAlmostEqual(score, 0.625)

    def test_markup_case_and_retracted_prefix_ignored(self):
        self.assertEqual(vr.dice("RETRACTED: FoxP3 Acts as a <i>Cotranscription</i> Factor",
                                 "FoxP3 acts as a cotranscription factor"), 1.0)

    def test_empty_or_stopword_titles_score_zero(self):
        self.assertEqual(vr.dice("", "PD-1 blockade"), 0.0)
        self.assertEqual(vr.dice("of the", "of the"), 0.0)

    def test_golden_dice_with_extra_words_on_both_sides(self):
        # A = {mutations, associated, acquired, resistance, pd, 1, blockade, melanoma} (8; with, to, in dropped)
        # B = {acquired, resistance, pd, 1, blockade, melanoma, mutations, jak1, jak2, b2m} (10; to, in, and dropped)
        # shared 7 -> Dice = 2*7/(8+10) = 14/18 = 0.7778 (recomputed by hand in review, November 2026)
        self.assertAlmostEqual(vr.dice("Mutations associated with acquired resistance to PD-1 blockade in melanoma",
                                       "Acquired resistance to PD-1 blockade in melanoma: mutations in JAK1, JAK2 and B2M"),
                               14 / 18)

    def test_duplicate_and_single_letter_words_collapse(self):
        # Word sets ignore repeats, and single letters are dropped: both titles reduce to {cells, tumours}.
        # Known limitation: titles that differ only in a single letter ('T cells' / 'B cells') score 1.0.
        self.assertEqual(vr.dice("T cells and T cells in tumours", "T cells in tumours"), 1.0)
        self.assertEqual(vr.dice("B cells in tumours", "T cells in tumours"), 1.0)


class AuthorsAndYears(unittest.TestCase):
    def test_surname_formats(self):
        for name in ("Ayers M", "M. Ayers", "Ayers, Mark", "Ayers M, Lunceford J", "Ayers MA et al."):
            self.assertEqual(vr.surname_key(name), "ayers", name)
        self.assertEqual(vr.surname_key("van der Burg SH"), "van der burg")
        self.assertEqual(vr.surname_key(""), "")

    def test_author_agreement(self):
        self.assertEqual(vr.compare_author("van der Burg SH", ["van der Burg"]), "agree")
        self.assertEqual(vr.compare_author("García-López J", ["Garcia-Lopez"]), "agree")
        self.assertEqual(vr.compare_author("Ribas A", ["Ayers", "Ayers M"]), "disagree")
        self.assertEqual(vr.compare_author("", ["Ayers"]), "not_compared")
        self.assertEqual(vr.compare_author("Ayers M", []), "not_compared")

    def test_first_author_from_list(self):
        self.assertEqual(vr.first_from_list("Ayers M, Lunceford J, Nebozhyn M"), "Ayers M")
        self.assertEqual(vr.first_from_list("Ayers, Mark; Lunceford, Jared"), "Ayers, Mark")

    def test_journal_abbreviations(self):
        self.assertTrue(vr.journals_agree("J Clin Invest", "Journal of Clinical Investigation"))
        self.assertTrue(vr.journals_agree("Nat Commun", "Nature communications"))
        self.assertTrue(vr.journals_agree("Proc Natl Acad Sci USA", "Proc Natl Acad Sci U S A"))
        self.assertFalse(vr.journals_agree("Nature", "Nature communications"))
        self.assertFalse(vr.journals_agree("Cell", "Cell Reports"))
        self.assertEqual(vr.compare_journal("", ["Nature"]), "not_compared")

    def test_year_checks(self):
        self.assertEqual(vr.compare_year(2017, [2016, 2017]), "exact")
        self.assertEqual(vr.compare_year(2016, [2017]), "off_by_one")
        self.assertEqual(vr.compare_year(2015, [2017, 2018]), "disagree")
        self.assertEqual(vr.compare_year(None, [2017]), "not_compared")
        self.assertEqual(vr.compare_year(2017, []), "not_compared")


class Notices(unittest.TestCase):
    def test_retraction_watch_record_in_crossref(self):
        rec = vr.crossref_record(CR_KURTULUS)
        self.assertEqual(len(rec["notices"]), 1)
        n = rec["notices"][0]
        self.assertEqual((n["kind"], n["source"], n["date"], n["notice"]),
                         ("retraction", "crossref/retraction-watch", "2024-07-01", "10.1172/jci183278"))
        self.assertEqual(vr.notice_status(rec["notices"], rec["notice_for"], vr.evidence_flags(rec, None)), "retracted")

    def test_article_replaced_by_its_notice_is_deduplicated(self):
        rec = vr.crossref_record(CR_IMMUNITY)
        self.assertEqual(sorted((n["source"], n["notice"]) for n in rec["notices"]),
                         [("crossref/publisher", "10.1016/j.immuni.2013.11.005"),
                          ("crossref/retraction-watch", "10.1016/j.immuni.2021.08.008")])
        self.assertEqual(rec["notice_for"], [])
        self.assertIn("title_says_retracted", vr.evidence_flags(rec, None))

    def test_cited_item_is_itself_a_notice(self):
        rec = vr.crossref_record(CR_NOTICE)
        flags = vr.evidence_flags(rec, None)
        self.assertEqual(rec["notice_for"][0]["target"], "10.5555/synthetic.article")
        self.assertIn("title_says_notice", flags)
        self.assertEqual(vr.notice_status(rec["notices"], rec["notice_for"], flags), "is_notice")

    def test_correction_and_concern_ranking(self):
        rec = vr.crossref_record(CR_CORRECTED)
        self.assertEqual(vr.notice_status(rec["notices"], rec["notice_for"], []), "corrected")
        concern = [{"kind": "concern"}, {"kind": "correction"}]
        self.assertEqual(vr.notice_status(concern, [], []), "expression_of_concern")
        self.assertEqual(vr.notice_status(concern + [{"kind": "retraction"}], [], []), "retracted")
        self.assertEqual(vr.notice_status([], [], []), "none_reported")

    def test_pubmed_retraction_not_double_counted(self):
        rec = vr.pubmed_record(PM_KURTULUS)
        self.assertEqual(len(rec["notices"]), 1)
        self.assertEqual(rec["notices"][0]["notice"], "PMID 38949028")
        self.assertEqual(rec["doi"], "10.1172/jci81187")

    def test_pubmed_notice_without_pmid_keeps_citation(self):
        # SYNTHETIC structure: PubMed sometimes lists a linked notice without a PMID.
        rec = vr.pubmed_record({"uid": "2", "title": "An example article.",
                                "references": [{"refsource": "Example J. 2014;1:e2.", "reftype": "Erratum in"}]})
        self.assertEqual(rec["notices"][0]["notice"], "Example J. 2014;1:e2.")
        self.assertEqual(vr.notice_status(rec["notices"], rec["notice_for"], []), "corrected")

    def test_pubmed_retraction_notice_record(self):
        rec = vr.pubmed_record({"uid": "1", "title": "Retraction.", "pubtype": ["Retraction of Publication"]})
        self.assertEqual(vr.notice_status(rec["notices"], rec["notice_for"], []), "is_notice")
        self.assertIsNone(vr.pubmed_record({"uid": "99999999", "error": "cannot get document summary"}))

    def test_current_pubmed_type_retraction_notice(self):
        # Live record PMID 38949028 (the Kurtulus retraction notice) carried the publication type
        # 'Retraction Notice' in November 2026 [NLM2026]; the type alone must mark the item as a notice.
        rec = vr.pubmed_record({"uid": "38949028", "title": "TIGIT predominantly regulates the immune response via regulatory T cells.",
                                "pubtype": ["Journal Article", "Retraction Notice"]})
        self.assertEqual(vr.notice_status(rec["notices"], rec["notice_for"], []), "is_notice")

    def test_notice_titles_versus_research_titles(self):
        for title in ("Correction to: An example article", "Retraction notice to: An example article",
                      "Expression of concern: An example article", "Correction for Smith et al., An example article"):
            self.assertTrue(vr.NOTICE_TITLE_RE.match(title), title)
        # SYNTHETIC research title: 'Correction for' followed by a method is not a notice.
        self.assertIsNone(vr.NOTICE_TITLE_RE.match("Correction for batch effects in single-cell RNA-seq data"))


class EvidenceFlags(unittest.TestCase):
    def test_review_and_trial_types(self):
        self.assertEqual(vr.evidence_flags(None, vr.pubmed_record(PM_CHEN)), ["review"])
        self.assertEqual(vr.evidence_flags(None, vr.pubmed_record(PM_AYERS)), ["clinical_trial"])

    def test_preprint_with_published_version(self):
        rec = vr.crossref_record(CR_PREPRINT)
        self.assertEqual(rec["published_versions"], ["10.2196/53164"])
        self.assertEqual(vr.evidence_flags(rec, None), ["preprint", "published_version_available"])

    def test_trial_protocol_is_not_trial_results(self):
        rec = vr.pubmed_record({"uid": "3", "title": "An example protocol.", "pubtype": ["Clinical Trial Protocol"]})
        self.assertEqual(vr.evidence_flags(None, rec), ["trial_protocol"])


class ClaimHelpers(unittest.TestCase):
    def test_numbers_skip_gene_names_and_years(self):
        self.assertEqual(vr.numbers_in("CD8, PD-1 and IL-2 in 2012"), [])
        self.assertEqual([v for _, v in vr.numbers_in("1,234 patients, 18% and 0.62")], [1234.0, 18.0, 0.62])

    def test_counts_in_year_range_are_kept_when_written_with_comma(self):
        self.assertEqual([v for _, v in vr.numbers_in("in 1,950 patients since 2012")], [1950.0])

    def test_both_ends_of_a_range_are_read(self):
        # Peters2025 abstract wording: '26-73% of cases' and '95% CI [3.06, 7.70]'; hyphen and en dash alike.
        self.assertEqual([v for _, v in vr.numbers_in("OR 4.85, 95% CI 3.06-7.70")], [4.85, 95.0, 3.06, 7.7])
        self.assertEqual([v for _, v in vr.numbers_in("in 26-73% of cases")], [26.0, 73.0])
        self.assertEqual([v for _, v in vr.numbers_in("in 26–73% of cases")], [26.0, 73.0])
        self.assertEqual(vr.claim_numbers_missing("overgeneralised in 73% of cases", "in 26-73% of cases"), [])

    def test_numbers_glued_to_letters_are_skipped(self):
        self.assertEqual(vr.numbers_in("TGF-β1, PD-1, IL-2, CD8, B16-F10, COVID-19, H-2Kb"), [])

    def test_claim_number_missing_from_abstract(self):
        claim = "Anti-PD-1 produced objective responses in 40% of patients with non-small-cell lung cancer."
        self.assertEqual(vr.claim_numbers_missing(claim, vr.clean_text(TOPALIAN_ABSTRACT)), ["40%"])
        self.assertEqual(vr.claim_numbers_missing("28% of melanoma patients responded", TOPALIAN_ABSTRACT), [])

    def test_model_cues(self):
        self.assertEqual(vr.model_cues("We compared melanoma growth in mice harboring distinct microbiota"), ["mouse"])
        self.assertEqual(vr.model_cues(TOPALIAN_ABSTRACT), ["human"])
        self.assertEqual(vr.model_cues(""), [])

    def test_murine_tumour_lines_and_strains_are_mouse_cues(self):
        for text in ("MC38 and CT26 tumours", "B16-F10 melanoma", "4T1 mammary carcinoma", "C57BL/6 hosts",
                     "GL261 glioma", "MC-38 tumours"):
            self.assertEqual(vr.model_cues(text), ["mouse"], text)
        self.assertEqual(vr.model_cues("CD4T1 and B2M"), [])


class Verdicts(unittest.TestCase):
    def setUp(self):
        self.ayers = vr.crossref_record(CR_AYERS)
        self.groom = vr.crossref_record(CR_GROOM)
        self.title = "IFN-gamma-related mRNA profile predicts clinical response to PD-1 blockade"

    def test_verified(self):
        r = vr.assess(ref(doi="10.1172/jci91190", title=self.title, year=2017, first_author="Ayers M"), self.ayers)
        self.assertEqual((r["verdict"], r["notice_status"], r["title_score"]), ("verified", "none_reported", 1.0))

    def test_wrong_first_author_is_metadata_error(self):
        r = vr.assess(ref(doi="10.1172/jci91190", title=self.title, year=2016, first_author="Ribas A"), self.ayers)
        self.assertEqual((r["verdict"], r["year_check"], r["author_check"]), ("metadata_errors", "off_by_one", "disagree"))
        self.assertTrue(any("Ribas" in i for i in r["issues"]))

    def test_doi_of_another_work_is_mismatch(self):
        r = vr.assess(ref(doi="10.1038/s41590-021-00999-x", year=2021, first_author="Lee K",
                          title="CXCL13+ exhausted CD8+ T cells predict anti-PD-1 response in hepatocellular carcinoma"), self.groom)
        self.assertEqual(r["verdict"], "mismatch")
        self.assertLess(r["title_score"], vr.TITLE_PARTIAL)
        # The record is the Groom editorial, not the cited work, so its notice status says nothing about the citation.
        self.assertEqual(r["notice_status"], "unknown")

    def test_title_thresholds_at_the_boundaries(self):
        # SYNTHETIC titles built by hand: 9 shared words of 10 + 10 gives 18/20 = 0.9 exactly (a tie with
        # TITLE_MATCH, counted as a match); 3 shared of 5 + 5 gives 0.6 (a tie with TITLE_PARTIAL, partial).
        a = "tumour immune cells express checkpoint ligands under interferon gamma stimulation"
        b = "tumour immune cells express checkpoint ligands under interferon gamma signalling"
        rec = dict(self.ayers, title=b, years=[2020], first_author="Smith")
        self.assertEqual(vr.assess(ref(doi="10.1/x", title=a, year=2020, first_author="Smith"), rec)["verdict"], "verified")
        rec = dict(rec, title="alpha beta gamma zeta theta")
        self.assertEqual(vr.assess(ref(doi="10.1/x", title="alpha beta gamma delta epsilon", year=2020), rec)["verdict"],
                         "check_title")
        rec = dict(rec, title="alpha beta zeta theta")      # 2 shared of 4 + 4 = 0.5
        self.assertEqual(vr.assess(ref(doi="10.1/x", title="alpha beta gamma delta", year=2020), rec)["verdict"], "mismatch")

    def test_doi_registered_outside_crossref_is_unconfirmed(self):
        # e.g. a DataCite DOI (arXiv 10.48550/...): Crossref answers 404, doi.org resolves it (checked November 2026).
        r = vr.assess(ref(doi="10.48550/arxiv.1706.03762"), doi_status="elsewhere")
        self.assertEqual((r["verdict"], r["notice_status"]), ("unconfirmed", "unknown"))
        self.assertTrue(any("not registered with Crossref" in i for i in r["issues"]))

    def test_mismatch_without_title_needs_two_disagreements(self):
        both = vr.assess(ref(doi="10.1038/s41590-021-00999-x", year=2015, first_author="Lee K"), self.groom)
        one = vr.assess(ref(doi="10.1038/s41590-021-00999-x", year=2021, first_author="Lee K"), self.groom)
        self.assertEqual((both["verdict"], one["verdict"]), ("mismatch", "metadata_errors"))

    def test_pmid_of_another_paper_without_title(self):
        # cited as Simoni et al., Nature 2018 but the PMID is Duhen et al., Nat Commun 2018: author and journal differ
        r = vr.assess(ref(pmid="30006565", year=2018, first_author="Simoni Y", journal="Nature"),
                      None, vr.pubmed_record(PM_DUHEN))
        self.assertEqual((r["verdict"], r["year_check"], r["journal_check"]), ("mismatch", "exact", "disagree"))

    def test_loosely_similar_search_hit_is_not_substituted(self):
        # Live case (Nov 2026): the invented 'Lee K 2021 Nat Immunol' title searched at Crossref scored
        # Dice 0.615 against a 2024 Cancers paper by Vahidian; year and author disagree, so it is rejected.
        hit = {"title": "CXCL13 Positive Cells Localization Predict Response to Anti-PD-1/PD-L1 in Pulmonary "
                        "Non-Small Cell Carcinoma", "years": [2024], "first_author": "Vahidian"}
        cited = ref(title="CXCL13+ exhausted CD8+ T cells predict anti-PD-1 response in hepatocellular carcinoma",
                    year=2021, first_author="Lee K")
        score = vr.dice(cited["title"], hit["title"])
        self.assertTrue(vr.TITLE_PARTIAL <= score < vr.TITLE_MATCH)
        self.assertFalse(vr.search_hit_acceptable(cited, score, hit))
        self.assertTrue(vr.search_hit_acceptable(dict(cited, year=2024, first_author="Vahidian F"), score, hit))
        self.assertFalse(vr.search_hit_acceptable(dict(cited, year=2024, first_author=""), score, hit))
        self.assertTrue(vr.search_hit_acceptable(cited, 0.95, hit))
        self.assertFalse(vr.search_hit_acceptable(cited, 0.5, hit))

    def test_identifier_only_is_unconfirmed(self):
        self.assertEqual(vr.assess(ref(doi="10.1172/jci91190"), self.ayers)["verdict"], "unconfirmed")

    def test_partial_title_needs_checking(self):
        r = vr.assess(ref(doi="10.1172/jci91190", title="IFN-gamma signature predicts PD-1 blockade response",
                          year=2017, first_author="Ayers"), self.ayers)
        self.assertEqual(r["verdict"], "check_title")

    def test_doi_and_pmid_of_different_works(self):
        r = vr.assess(ref(doi="10.1172/jci91190", pmid="28102259", title=self.title), self.ayers, vr.pubmed_record(PM_CHEN))
        self.assertEqual(r["verdict"], "identifier_conflict")

    def test_broken_doi_found_by_title(self):
        r = vr.assess(ref(doi="10.1172/jci9119", title=self.title, year=2017, first_author="Ayers"), self.ayers,
                      matched_by="title_search", doi_status="not_found")
        self.assertEqual(r["verdict"], "wrong_identifier")

    def test_nothing_found_versus_lookup_failure(self):
        self.assertEqual(vr.assess(ref(doi="10.1/x"), doi_status="not_found")["verdict"], "not_found")
        failed = vr.assess(ref(doi="10.1/x"), errors=["HTTP 503 from api.crossref.org"])
        self.assertEqual((failed["verdict"], failed["notice_status"]), ("lookup_failed", "unknown"))


class InputParsing(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def test_csv_aliases_bom_and_identifier_cleanup(self):
        path = write(self.dir, "refs.csv", "﻿ID,First Author,Year,Title,DOI,PMID\n"
                     "r1,Ayers M,2017,IFN profile,https://doi.org/10.1172/JCI91190.,PMID: 28650338\n")
        (r,) = vr.read_references(path)
        self.assertEqual((r["id"], r["doi"], r["pmid"], r["year"], r["first_author"]),
                         ("r1", "10.1172/jci91190", "28650338", 2017, "Ayers M"))

    def test_json_object_tsv_and_authors_list(self):
        path = write(self.dir, "refs.json", json.dumps({"references": [{"title": "X", "year": 2017,
                                                                        "authors": ["Ayers M", "Lunceford J"]}]}))
        self.assertEqual(vr.read_references(path)[0]["first_author"], "Ayers M")
        tsv = write(self.dir, "refs.tsv", "doi\ttitle\n10.1172/JCI91190\tIFN, profile\n")
        self.assertEqual(vr.read_references(tsv)[0]["title"], "IFN, profile")

    def test_bad_values_are_reported_not_fatal(self):
        path = write(self.dir, "refs.csv", "doi,title,year\nnot-a-doi,Some title,20177\n")
        (r,) = vr.read_references(path)
        self.assertEqual((r["doi"], r["year"]), ("", None))
        self.assertEqual(len(r["input_issues"]), 2)

    def test_common_citation_formatting_is_accepted(self):
        self.assertEqual(vr.normalise_doi("doi.org/10.1172/JCI91190"), ("10.1172/jci91190", None))
        self.assertEqual(vr.normalise_doi("10.1016/0895-4356(95)00048-8"), ("10.1016/0895-4356(95)00048-8", None))
        self.assertEqual(vr.normalise_pmid("PMID: 28650338."), ("28650338", None))
        for raw in ("2017", "2017a", "2017 Jun 26", "2017/06/26/"):
            self.assertEqual(vr.normalise_year(raw), (2017, None), raw)
        self.assertIsNone(vr.normalise_year("c2017")[0])

    def test_non_utf8_file_is_an_input_error(self):
        path = os.path.join(self.dir, "latin1.csv")
        with open(path, "wb") as handle:
            handle.write("doi,title\n10.1172/JCI91190,Lavou\xe9\n".encode("latin-1"))
        with self.assertRaises(vr.InputError):
            vr.read_references(path)

    def test_fatal_input_errors(self):
        cases = {"none.csv": "doi,title\n,\n", "cols.csv": "name,page\nA,1\n", "empty.json": "[]",
                 "broken.json": "[{"}
        for name, text in cases.items():
            with self.assertRaises(vr.InputError, msg=name):
                vr.read_references(write(self.dir, name, text))
        with self.assertRaises(vr.InputError):
            vr.read_references(os.path.join(self.dir, "missing.csv"))

    def test_cache_names(self):
        self.assertEqual(vr.cache_name("crossref-work", "10.1172/JCI91190"), "crossref-work--10.1172_jci91190.json")
        long_name = vr.cache_name("crossref-search", "word " * 60)
        self.assertLess(len(long_name), 120)
        self.assertNotEqual(long_name, vr.cache_name("crossref-search", "word " * 61))


class OfflineLookup(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name
        save(self.dir, "crossref-work", "10.1172/jci91190", {"message": CR_AYERS})
        save(self.dir, "pubmed-esearch", "10.1172/jci91190", esearch("28650338"))
        save(self.dir, "pubmed-esummary", "28650338", esummary(PM_AYERS))
        save(self.dir, "crossref-work", "10.1172/jci81187", {"message": CR_KURTULUS})
        save(self.dir, "pubmed-esearch", "10.1172/jci81187", esearch("26413872"))
        save(self.dir, "pubmed-esummary", "26413872", esummary(PM_KURTULUS))
        save(self.dir, "crossref-work", "10.1038/s41590-021-00999-x", {"message": CR_GROOM})
        save(self.dir, "pubmed-esearch", "10.1038/s41590-021-00999-x", esearch())
        # Title search for the invented 'Lee K' citation (its DOI belongs to Groom 2021): nothing similar exists.
        self.lee_title = "CXCL13+ exhausted CD8+ T cells predict anti-PD-1 response in hepatocellular carcinoma"
        save(self.dir, "crossref-search", f"{self.lee_title} Lee K 2021", {"message": {"items": [CR_GROOM, CR_AYERS]}})
        self.fetcher = vr.Fetcher(offline_dir=self.dir)

    def tearDown(self):
        self.tmp.cleanup()

    def test_doi_lookup_reaches_pubmed(self):
        r = vr.lookup(ref(doi="10.1172/jci91190", year=2017, first_author="Ayers M",
                          title="IFN-gamma-related mRNA profile predicts clinical response to PD-1 blockade"), self.fetcher)
        self.assertEqual((r["verdict"], r["record"]["pmid"]), ("verified", "28650338"))

    def test_retracted_paper_with_correct_metadata(self):
        r = vr.lookup(ref(doi="10.1172/jci81187", year=2015, first_author="Kurtulus S",
                          title="TIGIT predominantly regulates the immune response via regulatory T cells"), self.fetcher)
        self.assertEqual((r["verdict"], r["notice_status"]), ("verified", "retracted"))
        self.assertEqual({n["source"] for n in r["notices"]}, {"crossref/retraction-watch", "pubmed"})

    def test_unresolvable_doi_falls_back_to_title_search(self):
        title = "CXCL13+ exhausted CD8+ T cells predict anti-PD-1 response in hepatocellular carcinoma"
        r0 = ref(doi="10.1038/s41590-021-09999-x", title=title, year=2021, first_author="Lee K")
        save(self.dir, "crossref-work", r0["doi"], {"__status": 404})
        save(self.dir, "doi-handle", r0["doi"], {"__status": 404})
        save(self.dir, "crossref-search", f"{title} Lee K 2021", {"message": {"items": [CR_GROOM, CR_AYERS]}})
        r = vr.lookup(r0, self.fetcher)
        self.assertEqual((r["verdict"], r["notice_status"]), ("not_found", "unknown"))
        self.assertTrue(any("does not resolve" in i for i in r["issues"]))
        self.assertTrue(any("closest title found" in i for i in r["issues"]))

    def test_identifier_of_unrelated_work_and_title_found_nowhere(self):
        r = vr.lookup(ref(doi="10.1038/s41590-021-00999-x", title=self.lee_title, year=2021, first_author="Lee K"),
                      self.fetcher)
        self.assertEqual(r["verdict"], "mismatch")
        self.assertTrue(any("no record with the cited title" in i for i in r["issues"]))
        self.assertFalse(any("lookup error" in i for i in r["issues"]))

    def test_identifier_of_unrelated_work_but_cited_title_exists(self):
        title = "IFN-gamma-related mRNA profile predicts clinical response to PD-1 blockade"
        save(self.dir, "crossref-search", f"{title} Ayers M 2017", {"message": {"items": [CR_GROOM, CR_AYERS]}})
        r = vr.lookup(ref(doi="10.1038/s41590-021-00999-x", title=title, year=2017, first_author="Ayers M"), self.fetcher)
        self.assertEqual((r["verdict"], r["matched_by"]), ("wrong_identifier", "title_search"))
        self.assertEqual((r["record"]["doi"], r["record"]["pmid"]), ("10.1172/jci91190", "28650338"))
        self.assertTrue(any("belongs to a different work: 'Diversity in science" in i for i in r["issues"]))

    def test_missing_saved_response_is_a_lookup_failure(self):
        r = vr.lookup(ref(doi="10.1056/nejmoa1200690", title="Safety, activity"), self.fetcher)
        self.assertEqual(r["verdict"], "lookup_failed")
        self.assertTrue(any("crossref-work--10.1056_nejmoa1200690.json" in i for i in r["issues"]))

    def test_abstract_and_claim_numbers(self):
        save(self.dir, "crossref-work", "10.1056/nejmoa1200690", {"message": {
            "DOI": "10.1056/nejmoa1200690", "type": "journal-article", "issued": {"date-parts": [[2012]]},
            "title": ["Safety, Activity, and Immune Correlates of Anti–PD-1 Antibody in Cancer"],
            "author": [{"family": "Topalian", "sequence": "first"}]}})
        save(self.dir, "pubmed-esearch", "10.1056/nejmoa1200690", esearch())
        save(self.dir, "europepmc", "doi:10.1056/nejmoa1200690", EPMC_TOPALIAN)
        r = vr.lookup(ref(doi="10.1056/nejmoa1200690", year=2012, first_author="Topalian SL",
                          title="Safety, activity, and immune correlates of anti-PD-1 antibody in cancer",
                          claim="Objective responses in 40% of patients with non-small-cell lung cancer"),
                      self.fetcher, with_abstracts=True)
        self.assertEqual(r["verdict"], "verified")
        self.assertIn("Cumulative response rates", r["abstract"])
        self.assertEqual(r["claim_numbers_not_in_abstract"], ["40%"])
        self.assertEqual((r["pmcid"], r["model_cues"]), ("PMC3544539", ["human"]))


class CommandLine(unittest.TestCase):
    def setUp(self):
        OfflineLookup.setUp(self)   # same temporary folder and saved responses
        self.csv = write(self.dir, "refs.csv",
                         "id,first_author,year,title,journal,doi\n"
                         "1,Ribas A,2016,IFN-gamma-related mRNA profile predicts clinical response to PD-1 blockade,J Clin Invest,10.1172/JCI91190\n"
                         "2,Kurtulus S,2015,TIGIT predominantly regulates the immune response via regulatory T cells,J Clin Invest,10.1172/JCI81187\n"
                         "3,Lee K,2021,CXCL13+ exhausted CD8+ T cells predict anti-PD-1 response in hepatocellular carcinoma,Nat Immunol,10.1038/s41590-021-00999-x\n")

    def tearDown(self):
        self.tmp.cleanup()

    def run_cli(self, *args):
        return subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True, encoding="utf-8", timeout=60)

    def test_offline_json_report(self):
        proc = self.run_cli(self.csv, "--offline", self.dir)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        report = json.loads(proc.stdout)
        self.assertEqual([r["verdict"] for r in report["results"]], ["metadata_errors", "verified", "mismatch"])
        self.assertEqual(report["summary"]["notice_status"], {"retracted": 1, "none_reported": 1, "unknown": 1})
        self.assertEqual(report["mode"], "offline")

    def test_offline_tsv_report(self):
        out = os.path.join(self.dir, "report.tsv")
        proc = self.run_cli(self.csv, "--offline", self.dir, "--format", "tsv", "--out", out)
        self.assertEqual(proc.returncode, 0, proc.stderr)
        with open(out, encoding="utf-8") as handle:
            lines = handle.read().splitlines()
        self.assertEqual(lines[0].split("\t")[:3], ["id", "verdict", "notice_status"])
        self.assertEqual(len(lines), 4)
        self.assertIn("retraction (crossref/retraction-watch, 2024-07-01", lines[2])

    def test_input_errors_exit_2(self):
        bad = write(self.dir, "bad.csv", "name,page\nA,1\n")
        proc = self.run_cli(bad, "--offline", self.dir)
        self.assertEqual(proc.returncode, 2)
        self.assertIn("no doi, pmid or title column", proc.stderr)
        self.assertEqual(self.run_cli(self.csv, "--offline", self.dir, "--save", self.dir).returncode, 2)
        # Output and save paths are checked before any lookup, so a long online run is not lost at the end.
        proc = self.run_cli(self.csv, "--offline", self.dir, "--out", os.path.join(self.dir, "no-such-folder", "r.json"))
        self.assertEqual(proc.returncode, 2)
        self.assertIn("--out folder does not exist", proc.stderr)
        self.assertEqual(self.run_cli(self.csv, "--save", self.csv).returncode, 2)
        self.assertEqual(self.run_cli(self.csv, "--offline", self.dir, "--delay", "-1").returncode, 2)

    def test_help(self):
        proc = self.run_cli("--help")
        self.assertEqual(proc.returncode, 0)
        self.assertIn("--offline", proc.stdout)


if __name__ == "__main__":
    unittest.main()
