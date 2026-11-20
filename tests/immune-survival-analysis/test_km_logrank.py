"""Tests for skills/immune-survival-analysis/scripts/km_logrank.py (run: python3 -m unittest discover -s tests/immune-survival-analysis -v)."""
import importlib.util
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "..", "skills", "immune-survival-analysis", "scripts", "km_logrank.py")
spec = importlib.util.spec_from_file_location("km_logrank", SCRIPT)
km = importlib.util.module_from_spec(spec)
spec.loader.exec_module(km)

# The 'aml' data of the R survival package (time in weeks, status 1 = relapse), as used for the
# R-validated tests of the owner's tmelens package. Reference values below are R output:
#   survdiff(Surv(time, status) ~ x, data = aml); survfit(Surv(time, status) ~ 1, data = aml)
AML_TIME = [9, 13, 13, 18, 23, 28, 31, 34, 45, 48, 161, 5, 5, 8, 8, 12, 16, 23, 27, 30, 33, 43, 45]
AML_STATUS = [1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1]
AML_X = ["Maintained"] * 11 + ["Nonmaintained"] * 12
Z975 = 1.959963984540054


class KaplanMeierTests(unittest.TestCase):
    def test_hand_example(self):
        # times 1, 2, 2+, 3, 4+ (+ = censored)
        # t=1: n=5, d=1 -> S = 4/5 = 0.8;            Greenwood sum = 1/(5*4) = 0.05
        # t=2: n=4, d=1 -> S = 0.8 * 3/4 = 0.6;      sum += 1/(4*3) -> 0.133333
        # t=3: n=2, d=1 -> S = 0.6 * 1/2 = 0.3;      sum += 1/(2*1) -> 0.633333
        # t=4: n=1, d=0 -> S stays 0.3
        table = km.kaplan_meier([1, 2, 2, 3, 4], [1, 1, 0, 1, 0])
        self.assertEqual([r["time"] for r in table], [1, 2, 3, 4])
        self.assertEqual([r["n_risk"] for r in table], [5, 4, 2, 1])
        for got, want in zip([r["survival"] for r in table], [0.8, 0.6, 0.3, 0.3]):
            self.assertAlmostEqual(got, want, places=12)
        self.assertAlmostEqual(table[2]["se"], 0.3 * math.sqrt(0.05 + 1 / 12 + 0.5), places=12)
        # log interval at t=1: 0.8 * exp(-1.96 * sqrt(0.05)); upper 0.8 * exp(+...) = 1.24 is capped at 1
        self.assertAlmostEqual(table[0]["ci_lower"], 0.8 * math.exp(-Z975 * math.sqrt(0.05)), places=10)
        self.assertEqual(table[0]["ci_upper"], 1.0)

    def test_median_and_its_interval_hand_example(self):
        # median = first time with S <= 0.5 -> 3. Lower band at t=2: 0.6 * exp(-1.96 * sqrt(0.133333))
        # = 0.293 <= 0.5, at t=1 it is 0.516 > 0.5 -> lower limit 2. The upper band is capped at 1
        # at every time, so the upper limit is undefined (None).
        med = km.median_survival(km.kaplan_meier([1, 2, 2, 3, 4], [1, 1, 0, 1, 0]))
        self.assertEqual(med, {"median": 3, "ci_lower": 2, "ci_upper": None})

    def test_matches_r_survfit_on_aml(self):
        table = km.kaplan_meier(AML_TIME, AML_STATUS)
        r_time = [5, 8, 9, 12, 13, 16, 18, 23, 27, 28]
        r_surv = [0.913043478261, 0.826086956522, 0.782608695652, 0.739130434783, 0.695652173913,
                  0.695652173913, 0.645962732919, 0.546583850932, 0.496894409938, 0.496894409938]
        self.assertEqual([r["time"] for r in table][:10], r_time)
        for got, want in zip([r["survival"] for r in table][:10], r_surv):
            self.assertAlmostEqual(got, want, places=10)
        self.assertEqual(km.median_survival(table)["median"], 27)

    def test_everyone_fails_at_once(self):
        table = km.kaplan_meier([5, 5, 5], [1, 1, 1])
        self.assertEqual(table[0]["survival"], 0.0)
        self.assertIsNone(table[0]["se"])
        med = km.median_survival(table)
        self.assertEqual((med["median"], med["ci_lower"], med["ci_upper"]), (5, 5, None))

    def test_empty_input(self):
        self.assertEqual(km.kaplan_meier([], []), [])
        self.assertEqual(km.median_survival([]), {"median": None, "ci_lower": None, "ci_upper": None})

    def test_reverse_km_follow_up(self):
        # flip the indicator: censorings at 2 and 4 become events. t=2: n=4 at risk, 1 'event'
        # -> 0.75; t=4: n=1, 1 'event' -> 0. First value <= 0.5 is at t=4.
        self.assertEqual(km.reverse_km_median([1, 2, 2, 3, 4], [1, 1, 0, 1, 0]), 4)

    def test_survival_at_beyond_follow_up(self):
        table = km.kaplan_meier([1, 2, 2, 3, 4], [1, 1, 0, 1, 0])
        early = km.survival_at(table, 0.5)
        self.assertEqual((early["survival"], early["n_risk"], early["beyond_follow_up"]), (1.0, 5, False))
        mid = km.survival_at(table, 2.5)
        self.assertAlmostEqual(mid["survival"], 0.6)
        self.assertEqual(mid["n_risk"], 2)
        self.assertTrue(km.survival_at(table, 10)["beyond_follow_up"])


class LogRankTests(unittest.TestCase):
    def test_hand_two_groups(self):
        # A fails at 1 and 3, B at 2 and 4; all events.
        # t=1: n=4 (2A), E_A = 1/2, V = 1*3/3 * 1/2 * 1/2 = 1/4
        # t=2: n=3 (1A), E_A = 1/3, V += 1*2/2 * 1/3 * 2/3 = 2/9
        # t=3: n=2 (1A), E_A = 1/2, V += 1*1/1 * 1/2 * 1/2 = 1/4
        # t=4: n=1, no contribution.  O_A - E_A = 2 - 4/3 = 2/3; V = 13/18; chi2 = (4/9) / (13/18) = 8/13
        res = km.logrank_test([1, 3, 2, 4], [1, 1, 1, 1], ["A", "A", "B", "B"])
        self.assertAlmostEqual(res["chi2"], 8 / 13, places=12)
        self.assertEqual(res["df"], 1)
        self.assertAlmostEqual(res["p"], math.erfc(math.sqrt(4 / 13)), places=12)

    def test_matches_r_survdiff_on_aml(self):
        res = km.logrank_test(AML_TIME, AML_STATUS, AML_X)
        self.assertAlmostEqual(res["chi2"], 3.39638869898, places=9)
        self.assertEqual([g["observed"] for g in res["groups"]], [7, 11])
        self.assertAlmostEqual(res["groups"][0]["expected"], 10.6893359923007, places=9)
        self.assertAlmostEqual(res["groups"][1]["expected"], 7.31066400769928, places=9)

    def test_three_groups_and_order_invariance(self):
        times = [3, 5, 7, 2, 9, 11, 4, 6, 8, 12, 1, 10]
        events = [1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1]
        groups = ["a", "b", "c"] * 4
        one = km.logrank_test(times, events, groups)
        two = km.logrank_test(times, events, groups, order=["c", "a", "b"])
        self.assertEqual(one["df"], 2)
        self.assertAlmostEqual(one["chi2"], two["chi2"], places=10)
        self.assertAlmostEqual(sum(g["observed"] - g["expected"] for g in one["groups"]), 0.0, places=10)

    def test_undefined_cases(self):
        self.assertIn("only one group", km.logrank_test([1, 2], [1, 0], ["x", "x"])["reason"])
        self.assertIn("no events", km.logrank_test([1, 2], [0, 0], ["x", "y"])["reason"])

    def test_chi2_tail_closed_forms(self):
        self.assertAlmostEqual(km.chi2_sf(Z975 ** 2, 1), 0.05, places=12)
        self.assertAlmostEqual(km.chi2_sf(-2 * math.log(0.05), 2), 0.05, places=12)  # df 2: exp(-x/2)
        self.assertAlmostEqual(km.chi2_sf(7.814727903251178, 3), 0.05, places=9)    # qchisq(0.95, 3)
        self.assertEqual(km.chi2_sf(0.0, 4), 1.0)


class SplitTests(unittest.TestCase):
    def test_median_ties_go_low(self):
        labels, order, info = km.split_score([1, 2, 2, 3], "median")
        self.assertEqual(labels, ["low", "low", "low", "high"])
        self.assertEqual(info["cut_values"], [2])

    def test_tertiles(self):
        labels, order, _ = km.split_score([1, 2, 3, 4, 5, 6], "tertile")
        self.assertEqual(order, ["low", "mid", "high"])
        self.assertEqual(labels, ["low", "low", "mid", "mid", "high", "high"])


class CommandLineTests(unittest.TestCase):
    def run_cli(self, text, *args):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "cohort.tsv")
            with open(path, "w", encoding="utf-8") as handle:
                handle.write(text)
            return subprocess.run([sys.executable, SCRIPT, "--input", path, *args], capture_output=True, text=True)

    def test_smoke_groups(self):
        rows = ["time\tevent\tgroup"] + [f"{t}\t{e}\t{g}" for t, e, g in zip(AML_TIME, AML_STATUS, AML_X)] + ["\t1\tMaintained"]
        proc = self.run_cli("\n".join(rows) + "\n", "--time", "time", "--event", "event", "--group", "group", "--times", "12,24")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        report = json.loads(proc.stdout)
        self.assertAlmostEqual(report["logrank"]["chi2"], 3.39638869898, places=9)
        self.assertEqual(report["input"]["rows_dropped_missing"], 1)
        self.assertEqual(len(report["groups"][0]["survival_at"]), 2)

    def test_smoke_score_split(self):
        rows = ["t\te\tscore"] + [f"{t}\t{e}\t{i % 7}" for i, (t, e) in enumerate(zip(AML_TIME, AML_STATUS))]
        proc = self.run_cli("\n".join(rows) + "\n", "--time", "t", "--event", "e", "--score", "score", "--split", "median")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(json.loads(proc.stdout)["grouping"]["method"], "median")

    def test_bad_event_coding_exits_2(self):
        proc = self.run_cli("time\tstatus\n5\tDead\n", "--time", "time", "--event", "status")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("must be 1 (event) or 0", proc.stderr)

    def test_missing_column_exits_2(self):
        proc = self.run_cli("time\tstatus\n5\t1\n", "--time", "os", "--event", "status")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("column 'os' not found", proc.stderr)


class TiedEventTimeTests(unittest.TestCase):
    """Hand-derived cases with several events at one time (added in the independent review)."""

    T = [2, 2, 2, 3, 5, 5, 7, 8, 8, 10]
    E = [1, 1, 0, 1, 1, 0, 0, 1, 1, 0]

    def test_km_with_tied_events_by_hand(self):
        # t=2: n=10, d=2, c=1 -> S = 0.8;                    Greenwood sum = 2/(10*8) = 0.025
        # t=3: n=7,  d=1      -> S = 0.8*6/7 = 0.685714;     sum += 1/(7*6)
        # t=5: n=6,  d=1, c=1 -> S = 0.685714*5/6 = 0.571429; sum += 1/(6*5)
        # t=7: n=4,  c=1      -> unchanged
        # t=8: n=3,  d=2      -> S = 0.571429/3 = 0.190476;  sum += 2/(3*1)
        table = km.kaplan_meier(self.T, self.E)
        s = [0.8, 0.8 * 6 / 7, 0.8 * 6 / 7 * 5 / 6, 0.8 * 6 / 7 * 5 / 6, 0.8 * 6 / 7 * 5 / 6 / 3, 0.8 * 6 / 7 * 5 / 6 / 3]
        for got, want in zip([r["survival"] for r in table], s):
            self.assertAlmostEqual(got, want, places=12)
        self.assertEqual([r["n_risk"] for r in table], [10, 7, 6, 4, 3, 1])
        self.assertAlmostEqual(table[4]["se"], s[4] * math.sqrt(2 / 80 + 1 / 42 + 1 / 30 + 2 / 3), places=12)
        # log-log interval at t=3: exp(-exp(log(-log S) +/- z * sqrt(sum) / |log S|))
        half = Z975 * math.sqrt(2 / 80 + 1 / 42) / abs(math.log(s[1]))
        centre = math.log(-math.log(s[1]))
        loglog = km.kaplan_meier(self.T, self.E, conf_type="log-log")
        self.assertAlmostEqual(loglog[1]["ci_lower"], math.exp(-math.exp(centre + half)), places=12)
        self.assertAlmostEqual(loglog[1]["ci_upper"], math.exp(-math.exp(centre - half)), places=12)

    def test_median_and_follow_up_with_tied_events(self):
        # median: first S <= 0.5 is at t=8. Lower band: t=2 0.8*exp(-z*sqrt(0.025)) = 0.587 > 0.5;
        # t=3 0.686*exp(-z*sqrt(0.0488)) = 0.445 -> 3. The upper band is >= 1 up to t=8 -> None.
        self.assertEqual(km.median_survival(km.kaplan_meier(self.T, self.E)), {"median": 8, "ci_lower": 3, "ci_upper": None})
        # reverse KM: censorings at 2, 5, 7, 10 become events -> 0.9, 0.75, 0.5625, then 0 at t=10
        self.assertEqual(km.reverse_km_median(self.T, self.E), 10)

    def test_logrank_with_tied_event_times_by_hand(self):
        # A: events at 2, 2, 5, censored at 8. B: event at 2, censored at 3, events at 5 and 6, censored at 9.
        # t=2: n=9 (4 A), d=3 (2 A): E_A = 4/3, V = 3*6/8 * 4/9 * 5/9 = 5/9
        # t=5: n=5 (2 A), d=2 (1 A): E_A = 4/5, V = 2*3/4 * 2/5 * 3/5 = 9/25
        # t=6: n=3 (1 A), d=1 (0 A): E_A = 1/3, V = 1*2/2 * 1/3 * 2/3 = 2/9
        # O_A - E_A = 3 - 37/15 = 8/15; V = 256/225; chi2 = (64/225) / (256/225) = 1/4
        res = km.logrank_test([2, 2, 5, 8, 2, 3, 5, 6, 9], [1, 1, 1, 0, 1, 0, 1, 1, 0], ["A"] * 4 + ["B"] * 5)
        self.assertAlmostEqual(res["chi2"], 0.25, places=12)
        self.assertAlmostEqual(res["groups"][0]["expected"], 37 / 15, places=12)
        self.assertAlmostEqual(res["p"], math.erfc(math.sqrt(0.125)), places=12)


class InputValidationTests(unittest.TestCase):
    run_cli = CommandLineTests.run_cli

    def test_duplicate_patient_id_exits_2(self):
        proc = self.run_cli("pid\ttime\tevent\nP1\t5\t1\nP2\t6\t0\nP1\t7\t1\n", "--time", "time", "--event", "event", "--id", "pid")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("'P1' (2 rows)", proc.stderr)

    def test_unique_patient_id_runs(self):
        proc = self.run_cli("pid\ttime\tevent\nP1\t5\t1\nP2\t6\t0\nP3\t7\t1\n", "--time", "time", "--event", "event", "--id", "pid")
        self.assertEqual(proc.returncode, 0, proc.stderr)
        self.assertEqual(json.loads(proc.stdout)["input"]["id_column_checked"], "pid")

    def test_negative_report_time_exits_2(self):
        proc = self.run_cli("time\tevent\n5\t1\n6\t0\n", "--time", "time", "--event", "event", "--times=-5")
        self.assertEqual(proc.returncode, 2)
        self.assertIn("--times must be finite and >= 0", proc.stderr)

    def test_unwritable_output_exits_2(self):
        with tempfile.TemporaryDirectory() as tmp:
            missing = os.path.join(tmp, "no_such_dir", "out.json")
            proc = self.run_cli("time\tevent\n5\t1\n6\t0\n", "--time", "time", "--event", "event", "--out", missing)
        self.assertEqual(proc.returncode, 2)
        self.assertIn("cannot write", proc.stderr)


if __name__ == "__main__":
    unittest.main()
