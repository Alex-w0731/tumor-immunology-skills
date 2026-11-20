#!/usr/bin/env python3
"""Validate every skill in this repository. Standard library only.

What is checked (CONTRIBUTING.md explains why each rule exists):

  frontmatter   name / description / license / metadata, length and charset rules
  body          line budget, relative links resolve, no orphan resource files
  evidence      references/sources.md exists, every [Key2020] used is defined there,
                every source entry carries a DOI, PMID or URL
  scripts       standard-library-only imports, argparse entry point, `--help` runs
  tests         every script has tests/<skill>/test_<script>.py
  evals         `claude plugin eval` case layout (prompt.md, case.yaml, graders/)
  hygiene       no local paths, e-mail addresses, TODO markers or tabs in shipped files
  manifests     .claude-plugin/plugin.json and marketplace.json

usage:
  python3 tools/validate_skills.py                 # everything
  python3 tools/validate_skills.py --skill NAME    # one skill
  python3 tools/validate_skills.py --no-run        # skip running `script --help`
  python3 tools/validate_skills.py --json          # machine-readable report

Exit code 0 when there are no errors. Warnings never fail the run.
"""
from __future__ import annotations

import argparse
import ast
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SKILLS_DIR = ROOT / "skills"
TESTS_DIR = ROOT / "tests"

NAME_RE = re.compile(r"^[a-z0-9]+(-[a-z0-9]+)*$")
MAX_NAME = 64
MAX_DESCRIPTION_SPEC = 1024      # open Agent Skills specification
MAX_DESCRIPTION_HOUSE = 900      # this repository's own budget
MAX_BODY_LINES = 500             # hard cap
WARN_BODY_LINES = 350
MIN_EVAL_CASES = 3
RESERVED_CASE_NAMES = {"results", "mocks", "node_modules", "files", "true", "false", "null"}
PROMPT_SETTING_KEYS = {"tags", "max_turns", "timeout_seconds", "allowed_tools"}
REQUIRED_METADATA = ("author", "version", "domain", "summary", "summary-zh")

LINK_RE = re.compile(r"(?<!\!)\[[^\]\n]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
IMAGE_RE = re.compile(r"\!\[[^\]\n]*\]\(([^)\s]+)\)")
SOURCE_KEY_RE = re.compile(r"\[([A-Z][A-Za-z'’\-]*(?:19|20)\d{2}[a-z]?)\]")
DOI_RE = re.compile(r"\b10\.\d{4,9}/(?:[^\s\"<>()\[\]]|\([^\s\"<>()]*\))+")
PMID_RE = re.compile(r"\bPMID:?\s*(\d{6,9})\b")
URL_RE = re.compile(r"https?://[^\s)>\]]+")
EMAIL_RE = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9-]+\.[A-Za-z]{2,}")
LOCAL_PATH_RE = re.compile(r"(/Users/[A-Za-z0-9._-]+/|/private/tmp/|/home/[a-z0-9._-]+/|[A-Z]:\\\\Users\\\\)")
MARKER_RE = re.compile(r"\b(TODO|FIXME|TBD|XXX|lorem ipsum)\b", re.IGNORECASE)

FALLBACK_STDLIB = {
    "__future__", "abc", "argparse", "array", "ast", "base64", "bisect", "collections", "contextlib",
    "copy", "csv", "dataclasses", "datetime", "decimal", "difflib", "enum", "fractions", "functools",
    "gzip", "hashlib", "heapq", "io", "itertools", "json", "math", "numbers", "operator", "os",
    "pathlib", "random", "re", "shutil", "statistics", "string", "subprocess", "sys", "tempfile",
    "textwrap", "typing", "unicodedata", "unittest", "warnings", "zipfile",
}
STDLIB = set(getattr(sys, "stdlib_module_names", ())) | FALLBACK_STDLIB


class Report:
    def __init__(self) -> None:
        self.items: list[dict] = []

    def add(self, level: str, where: str, message: str) -> None:
        self.items.append({"level": level, "where": where, "message": message})

    def error(self, where: str, message: str) -> None:
        self.add("error", where, message)

    def warn(self, where: str, message: str) -> None:
        self.add("warning", where, message)

    @property
    def errors(self) -> list[dict]:
        return [i for i in self.items if i["level"] == "error"]

    @property
    def warnings(self) -> list[dict]:
        return [i for i in self.items if i["level"] == "warning"]


# --------------------------------------------------------------------------- YAML subset

def _scalar(value: str):
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in "\"'":
        return value[1:-1]
    if value.startswith("[") and value.endswith("]"):
        inner = value[1:-1].strip()
        return [_scalar(v) for v in inner.split(",")] if inner else []
    return value


def parse_simple_yaml(text: str) -> dict:
    """Parse the small YAML subset used in frontmatter and case.yaml.

    Supports `key: value`, quoted scalars, inline lists, block scalars (`|`, `>` with
    optional `-`/`+`), and one level of nested mappings. Raises ValueError otherwise.
    """
    result: dict = {}
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if not line.strip() or line.lstrip().startswith("#"):
            i += 1
            continue
        if line[0] in " \t":
            raise ValueError(f"unexpected indentation at line {i + 1}: {line!r}")
        if ":" not in line:
            raise ValueError(f"expected 'key: value' at line {i + 1}: {line!r}")
        key, _, rest = line.partition(":")
        key, rest = key.strip(), rest.strip()
        if key in result:
            raise ValueError(f"duplicate key {key!r}")
        block: list[str] = []
        j = i + 1
        while j < len(lines) and (not lines[j].strip() or lines[j][0] in " \t"):
            block.append(lines[j])
            j += 1
        if rest[:1] in ("|", ">") and rest.rstrip("+-") in ("|", ">"):
            kept = [b for b in block]
            indents = [len(b) - len(b.lstrip()) for b in kept if b.strip()]
            pad = min(indents) if indents else 0
            body = [b[pad:] if b.strip() else "" for b in kept]
            while body and not body[-1]:
                body.pop()
            result[key] = "\n".join(body) if rest[0] == "|" else " ".join(s.strip() for s in body if s.strip())
        elif rest == "":
            nested: dict = {}
            items: list = []
            for b in block:
                if not b.strip() or b.lstrip().startswith("#"):
                    continue
                s = b.strip()
                if s.startswith("- "):
                    items.append(_scalar(s[2:]))
                elif ":" in s:
                    k, _, v = s.partition(":")
                    nested[k.strip()] = _scalar(v)
                else:
                    raise ValueError(f"cannot parse nested line under {key!r}: {b!r}")
            result[key] = items if items and not nested else nested
        else:
            if any(b.strip() for b in block):
                # plain multi-line scalar continuation
                rest = " ".join([rest] + [b.strip() for b in block if b.strip()])
            result[key] = _scalar(rest)
        i = j
    return result


def split_frontmatter(text: str) -> tuple[str | None, str]:
    if not text.startswith("---\n"):
        return None, text
    end = text.find("\n---", 4)
    if end == -1:
        return None, text
    after = text.find("\n", end + 1)
    return text[4:end], (text[after + 1:] if after != -1 else "")


# --------------------------------------------------------------------------- checks

def rel(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def check_links(md: Path, text: str, rep: Report) -> None:
    for regex in (LINK_RE, IMAGE_RE):
        for target in regex.findall(text):
            if re.match(r"^(https?:|mailto:|#)", target):
                continue
            path = (md.parent / target.split("#", 1)[0]).resolve()
            if not path.exists():
                rep.error(rel(md), f"broken relative link: {target}")


def strip_code(text: str) -> str:
    text = re.sub(r"```.*?```", "", text, flags=re.DOTALL)
    return re.sub(r"`[^`\n]*`", "", text)


def check_hygiene(path: Path, text: str, rep: Report) -> None:
    where = rel(path)
    if LOCAL_PATH_RE.search(text):
        rep.error(where, "contains a machine-specific absolute path")
    for m in EMAIL_RE.finditer(text):
        if not m.group(0).endswith(("example.com", "example.org")):
            rep.error(where, f"contains an e-mail address ({m.group(0).split('@')[0][:2]}…@…)")
            break
    if path.suffix == ".md" and MARKER_RE.search(strip_code(text)):
        rep.warn(where, "contains a TODO / FIXME / TBD marker")
    if text and not text.endswith("\n"):
        rep.warn(where, "does not end with a newline")


def check_plain_scalars(front: str, where: str, rep: Report) -> None:
    """Flag unquoted YAML values that a real YAML parser rejects or silently truncates.

    The light parser above accepts them, but Claude Code and other agents use full YAML, where
    an unquoted value containing ': ' is an error and ' #' starts a comment.
    """
    for n, line in enumerate(front.splitlines(), 1):
        m = re.match(r"^(\s*)([A-Za-z0-9_-]+):\s+(.*)$", line)
        if not m:
            continue
        value = m.group(3).rstrip()
        if not value or value[0] in "\"'" or value in ("|", ">", "|-", ">-", "|+", ">+"):
            continue
        if value.startswith("[") and value.endswith("]"):
            continue
        if value[0] in "&*!%@`{|>":
            rep.error(where, f"frontmatter line {n}: unquoted value starts with a YAML indicator character; wrap it in double quotes")
        elif re.search(r":\s|\s#", value) or value.endswith(":"):
            rep.error(where, f"frontmatter line {n}: unquoted value contains ': ' or ' #', which YAML rejects or truncates; wrap it in double quotes or use a '>-' block")


def check_frontmatter(skill: Path, rep: Report) -> tuple[dict, str]:
    md = skill / "SKILL.md"
    where = rel(md)
    text = md.read_text(encoding="utf-8")
    front, body = split_frontmatter(text)
    if front is None:
        rep.error(where, "missing YAML frontmatter delimited by '---' lines")
        return {}, text
    check_plain_scalars(front, where, rep)
    try:
        meta = parse_simple_yaml(front)
    except ValueError as exc:
        rep.error(where, f"frontmatter is not parseable: {exc}")
        return {}, body
    name = meta.get("name")
    if not isinstance(name, str) or not name:
        rep.error(where, "frontmatter needs a 'name'")
    else:
        if name != skill.name:
            rep.error(where, f"name {name!r} must equal the folder name {skill.name!r}")
        if not NAME_RE.match(name) or len(name) > MAX_NAME:
            rep.error(where, "name must be lowercase letters, digits and single hyphens, at most 64 characters")
        if re.search(r"(anthropic|claude)", name):
            rep.error(where, "name must not contain the reserved words 'anthropic' or 'claude'")
    desc = meta.get("description")
    if not isinstance(desc, str) or not desc.strip():
        rep.error(where, "frontmatter needs a non-empty 'description'")
    else:
        if len(desc) > MAX_DESCRIPTION_SPEC:
            rep.error(where, f"description is {len(desc)} characters; the specification allows {MAX_DESCRIPTION_SPEC}")
        elif len(desc) > MAX_DESCRIPTION_HOUSE:
            rep.warn(where, f"description is {len(desc)} characters; house budget is {MAX_DESCRIPTION_HOUSE}")
        if "<" in desc or ">" in desc:
            rep.error(where, "description must not contain angle brackets")
        if re.search(r"\b(I|we|you)\b", desc[:200]):
            rep.warn(where, "description should be written in the third person")
        if not re.search(r"\buse (this skill |it )?(when|whenever|for|to)\b", desc, re.IGNORECASE):
            rep.warn(where, "description should say when to use the skill ('Use when ...')")
    if meta.get("license") != "MIT":
        rep.error(where, "frontmatter 'license' must be MIT")
    md_meta = meta.get("metadata")
    if not isinstance(md_meta, dict):
        rep.error(where, "frontmatter needs a 'metadata' mapping")
    else:
        for key in REQUIRED_METADATA:
            if not str(md_meta.get(key, "")).strip():
                rep.error(where, f"metadata.{key} is required")
        for key in ("summary", "summary-zh"):
            if len(str(md_meta.get(key, ""))) > 160:
                rep.warn(where, f"metadata.{key} is longer than 160 characters")
        if not re.match(r"^\d+\.\d+\.\d+$", str(md_meta.get("version", ""))):
            rep.error(where, "metadata.version must look like 1.2.3 (quote it in YAML)")
    unknown = set(meta) - {"name", "description", "license", "metadata", "compatibility", "allowed-tools"}
    if unknown:
        rep.warn(where, f"unexpected frontmatter keys: {', '.join(sorted(unknown))}")
    return meta, body


def check_body(skill: Path, body: str, rep: Report) -> None:
    md = skill / "SKILL.md"
    where = rel(md)
    n = len(body.strip("\n").splitlines())
    if n > MAX_BODY_LINES:
        rep.error(where, f"body is {n} lines; the cap is {MAX_BODY_LINES}")
    elif n > WARN_BODY_LINES:
        rep.warn(where, f"body is {n} lines; aim for {WARN_BODY_LINES} or fewer and move detail to references/")
    if not re.search(r"^# \S", body, re.MULTILINE):
        rep.error(where, "body needs a level-one title")
    for heading in ("Scope", "Workflow", "Pitfalls", "Scripts", "Report", "References"):
        if not re.search(rf"^## .*{heading}", body, re.MULTILINE | re.IGNORECASE):
            rep.warn(where, f"expected a '## {heading}' section")
    check_links(md, body, rep)
    for sub in ("references", "scripts"):
        folder = skill / sub
        if not folder.is_dir():
            rep.error(rel(skill), f"missing {sub}/ folder")
            continue
        for f in sorted(folder.iterdir()):
            if f.name.startswith(".") or f.name == "__pycache__":
                continue
            if f.is_dir():
                rep.warn(rel(f), "nested folders are discouraged; keep resources one level deep")
                continue
            if f"{sub}/{f.name}" not in body:
                rep.error(where, f"{sub}/{f.name} is never mentioned in SKILL.md (orphan resource)")
    extras = {p.name for p in skill.iterdir()} - {"SKILL.md", "references", "scripts", "evals", "assets", ".DS_Store"}
    if extras:
        rep.warn(rel(skill), f"unexpected entries in the skill folder: {', '.join(sorted(extras))}")


def check_evidence(skill: Path, body: str, rep: Report) -> None:
    sources = skill / "references" / "sources.md"
    if not sources.is_file():
        rep.error(rel(skill), "references/sources.md is required (bibliography for every claim)")
        return
    src_text = sources.read_text(encoding="utf-8")
    defined: dict[str, str] = {}
    current = None
    for line in src_text.splitlines():
        m = re.match(r"^\s*(?:[-*]\s+)?(?:\*\*)?\[([A-Z][A-Za-z'’\-]*(?:19|20)\d{2}[a-z]?)\]", line)
        if m:
            current = m.group(1)
            defined[current] = line
        elif current and line.strip() and not line.startswith("#"):
            defined[current] += " " + line
        elif not line.strip() or line.startswith("#"):
            current = None
    if not defined:
        rep.error(rel(sources), "no source entries found; start each entry with its key, for example '- **[Ayers2017]** ...'")
    for key, entry in defined.items():
        if not (DOI_RE.search(entry) or PMID_RE.search(entry) or URL_RE.search(entry)):
            rep.error(rel(sources), f"source [{key}] has no DOI, PMID or URL")
    used: dict[str, str] = {}
    texts = {skill / "SKILL.md": body}
    for f in sorted((skill / "references").glob("*.md")):
        if f.name != "sources.md":
            texts[f] = f.read_text(encoding="utf-8")
    for path, text in texts.items():
        for key in SOURCE_KEY_RE.findall(text):
            used.setdefault(key, rel(path))
    for key, where in sorted(used.items()):
        if key not in defined:
            rep.error(where, f"source key [{key}] is not defined in references/sources.md")
    unused = sorted(set(defined) - set(used))
    if unused:
        rep.warn(rel(sources), f"sources never cited in the skill: {', '.join(unused)}")
    if not used:
        rep.warn(rel(skill), "SKILL.md and references cite no source keys")
    for f in sorted((skill / "references").glob("*.md")):
        text = f.read_text(encoding="utf-8")
        check_links(f, text, rep)
        check_hygiene(f, text, rep)
        n = len(text.splitlines())
        if n > 100 and not re.search(r"^## (Contents|Table of contents)", text, re.MULTILINE | re.IGNORECASE):
            rep.warn(rel(f), f"{n} lines without a '## Contents' section near the top")


def check_scripts(skill: Path, rep: Report, run: bool, meta: dict | None = None) -> None:
    scripts = sorted((skill / "scripts").glob("*.py")) if (skill / "scripts").is_dir() else []
    if not scripts:
        rep.error(rel(skill), "at least one script under scripts/ is required")
    for script in scripts:
        where = rel(script)
        text = script.read_text(encoding="utf-8")
        check_hygiene(script, text, rep)
        if not text.startswith("#!/usr/bin/env python3"):
            rep.error(where, "must start with '#!/usr/bin/env python3'")
        try:
            tree = ast.parse(text)
        except SyntaxError as exc:
            rep.error(where, f"syntax error: {exc}")
            continue
        if not ast.get_docstring(tree):
            rep.error(where, "needs a module docstring with usage examples")
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(a.name.split(".")[0] for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                imported.add(node.module.split(".")[0])
        third_party = sorted(m for m in imported if m not in STDLIB)
        if third_party:
            rep.error(where, f"imports outside the standard library: {', '.join(third_party)}")
        if "argparse" not in imported:
            rep.error(where, "must expose an argparse command line")
        if "__main__" not in text:
            rep.error(where, "needs an `if __name__ == \"__main__\":` entry point")
        # Network use is allowed only when the skill declares it in the `compatibility` field.
        declares_network = "network" in str((meta or {}).get("compatibility", "")).lower()
        for banned in ("urllib", "socket", "http", "requests"):
            if banned in imported and not declares_network:
                rep.error(where, f"imports {banned} but the frontmatter `compatibility` field does not declare network access")
        test = TESTS_DIR / skill.name / f"test_{script.stem}.py"
        if not test.is_file():
            rep.error(where, f"missing unit tests: {rel(test)}")
        if run:
            try:
                proc = subprocess.run([sys.executable, str(script), "--help"], capture_output=True, text=True, timeout=30)
                if proc.returncode != 0 or "usage" not in proc.stdout.lower():
                    rep.error(where, f"`--help` failed (exit {proc.returncode}): {proc.stderr.strip()[:200]}")
            except subprocess.TimeoutExpired:
                rep.error(where, "`--help` timed out")


def check_evals(skill: Path, rep: Report) -> None:
    evals = skill / "evals"
    if not evals.is_dir():
        rep.error(rel(skill), "missing evals/ folder")
        return
    cases = [d for d in sorted(evals.iterdir()) if d.is_dir() and d.name != "results"]
    if len(cases) < MIN_EVAL_CASES:
        rep.error(rel(evals), f"needs at least {MIN_EVAL_CASES} eval cases, found {len(cases)}")
    for case in cases:
        where = rel(case)
        if not NAME_RE.match(case.name) or not case.name[0].isalpha() or case.name in RESERVED_CASE_NAMES:
            rep.error(where, "invalid case name (lowercase letters, digits, hyphens; starts with a letter; not reserved)")
        prompt = case / "prompt.md"
        if not prompt.is_file():
            rep.error(where, "missing prompt.md")
        else:
            text = prompt.read_text(encoding="utf-8")
            check_hygiene(prompt, text, rep)
            front, body = split_frontmatter(text)
            if front is None:
                rep.error(rel(prompt), "prompt.md must start with a settings block")
            else:
                try:
                    settings = parse_simple_yaml(front)
                    extra = set(settings) - PROMPT_SETTING_KEYS
                    if extra:
                        rep.error(rel(prompt), f"unknown settings keys: {', '.join(sorted(extra))}")
                except ValueError as exc:
                    rep.error(rel(prompt), f"settings block not parseable: {exc}")
                if len(body.strip()) < 40:
                    rep.error(rel(prompt), "the prompt body is missing or too short to be realistic")
            if (case / "files").is_dir():
                if not (case / "place-files.sh").is_file():
                    rep.error(where, "files/ present but place-files.sh is missing")
                if "Input files:" not in body if front is not None else True:
                    rep.error(rel(prompt), "a case with files/ must end its prompt with an 'Input files: ...' line")
        case_yaml = case / "case.yaml"
        if not case_yaml.is_file():
            rep.error(where, "missing case.yaml")
        else:
            raw = case_yaml.read_text(encoding="utf-8")
            if not re.search(r'^schema_version:\s*"1\.1"\s*$', raw, re.MULTILINE):
                rep.error(rel(case_yaml), 'schema_version must be the quoted string "1.1"')
            m = re.search(r"^name:\s*(\S+)\s*$", raw, re.MULTILINE)
            if not m or m.group(1).strip("\"'") != case.name:
                rep.error(rel(case_yaml), "name must equal the case folder name")
            if not re.search(r"^expected_outcome:\s*\|\s*\n\s+\S", raw, re.MULTILINE):
                rep.error(rel(case_yaml), "expected_outcome must be a non-empty '|' block")
            if (case / "files").is_dir() and "scaffold_script: place-files.sh" not in raw:
                rep.error(rel(case_yaml), "a case with files/ must declare context.scaffold_script: place-files.sh")
        graders = sorted((case / "graders").glob("*.md")) if (case / "graders").is_dir() else []
        if not graders:
            rep.error(where, "needs at least one grader in graders/")
        elif len(graders) < 3:
            rep.warn(where, f"only {len(graders)} grader(s); three to six make a case informative")
        for g in graders:
            text = g.read_text(encoding="utf-8")
            front, body = split_frontmatter(text)
            if not re.match(r"^\d{2}-[a-z0-9-]+\.md$", g.name):
                rep.error(rel(g), "grader files are named NN-few-words.md")
            if front is None or not all(k in front for k in ("type:", "focus:", "weight:")):
                rep.error(rel(g), "grader needs the header block with type, focus and weight")
            elif len(body.strip()) < 15:
                rep.error(rel(g), "grader statement is missing")


def check_manifests(skill_names: list[str], rep: Report) -> None:
    plugin = ROOT / ".claude-plugin" / "plugin.json"
    market = ROOT / ".claude-plugin" / "marketplace.json"
    for path, required in ((plugin, ("name", "description", "version")), (market, ("name", "owner", "plugins"))):
        if not path.is_file():
            rep.error(rel(path), "file is missing")
            continue
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            rep.error(rel(path), f"invalid JSON: {exc}")
            continue
        for key in required:
            if key not in data:
                rep.error(rel(path), f"missing required key {key!r}")
        if path == market:
            plugins = data.get("plugins") or []
            if not isinstance(plugins, list) or not plugins:
                rep.error(rel(path), "plugins must be a non-empty list")
            for p in plugins if isinstance(plugins, list) else []:
                if not isinstance(p, dict) or "name" not in p or "source" not in p:
                    rep.error(rel(path), "each plugin entry needs 'name' and 'source'")
                    continue
                listed = p.get("skills")
                if isinstance(listed, list):
                    listed_names = sorted(Path(s).name for s in listed)
                    if listed_names != sorted(skill_names):
                        rep.error(rel(path), "the 'skills' list is out of sync with the skills/ folder")


def validate(only: str | None, run: bool) -> Report:
    rep = Report()
    if not SKILLS_DIR.is_dir():
        rep.error("skills/", "folder is missing")
        return rep
    skills = [d for d in sorted(SKILLS_DIR.iterdir()) if d.is_dir() and not d.name.startswith(".")]
    names = [s.name for s in skills]
    for skill in skills:
        if only and skill.name != only:
            continue
        if not (skill / "SKILL.md").is_file():
            rep.error(rel(skill), "missing SKILL.md")
            continue
        meta, body = check_frontmatter(skill, rep)
        check_hygiene(skill / "SKILL.md", (skill / "SKILL.md").read_text(encoding="utf-8"), rep)
        check_body(skill, body, rep)
        check_evidence(skill, body, rep)
        check_scripts(skill, rep, run, meta)
        check_evals(skill, rep)
    if not only:
        check_manifests(names, rep)
        if not names:
            rep.error("skills/", "no skills found")
    return rep


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate the skills in this repository.")
    parser.add_argument("--skill", help="validate a single skill by folder name")
    parser.add_argument("--no-run", action="store_true", help="do not execute `script --help`")
    parser.add_argument("--json", action="store_true", help="print a JSON report")
    args = parser.parse_args(argv)
    rep = validate(args.skill, run=not args.no_run)
    if args.json:
        print(json.dumps({"errors": rep.errors, "warnings": rep.warnings}, indent=2, ensure_ascii=False))
    else:
        for item in rep.items:
            print(f"{item['level'].upper():7} {item['where']}: {item['message']}")
        print(f"\n{len(rep.errors)} error(s), {len(rep.warnings)} warning(s)")
    return 1 if rep.errors else 0


if __name__ == "__main__":
    sys.exit(main())
