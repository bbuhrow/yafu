"""
yafu_schema.py -- extract yafu's command-line option table from cmdOptions.c.

yafu's options live in four parallel C arrays in top/cmdParser/cmdOptions.c:
    OptionArray[]       option names            ("B1ecm", "threads", ...)
    OptionHelp[]        help strings            ("(Integer < 64-bit): B1 bound ...")
    needsArg[]          0 = flag, 1 = required arg, 2 = optional arg
    LongOptionAliases[] (currently all empty)

Rather than hand-copying 149 entries into the GUI, this module parses those
arrays and writes yafu_options.json.  Re-run it whenever options are added:

    python yafu_schema.py path/to/top/cmdParser/cmdOptions.c  [-o yafu_options.json]

The GUI loads the JSON; its curated layout (yafu_gui.py: MODE_GROUPS) only
references option *names*, so help text and arg-ness always come from here.
"""
import json
import re
from pathlib import Path

_STR = re.compile(r'"((?:[^"\\]|\\.)*)"')


def _array_body(src: str, name: str) -> str:
    m = re.search(name + r'\s*\[[^\]]*\](?:\s*\[[^\]]*\])?\s*=\s*\{', src)
    if not m:
        raise ValueError(f"array {name} not found")
    i, depth = m.end(), 1
    j = i
    while depth:
        c = src[j]
        if c == '"':                      # skip string literal
            j += 1
            while src[j] != '"':
                j += 2 if src[j] == '\\' else 1
        elif c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
        j += 1
    return src[i:j - 1]


def _strip_comments(s: str) -> str:
    out, i, n = [], 0, len(s)
    while i < n:
        if s[i] == '"':
            j = i + 1
            while s[j] != '"':
                j += 2 if s[j] == '\\' else 1
            out.append(s[i:j + 1]); i = j + 1
        elif s.startswith('//', i):
            i = s.find('\n', i); i = n if i < 0 else i
        elif s.startswith('/*', i):
            i = s.find('*/', i) + 2
        else:
            out.append(s[i]); i += 1
    return ''.join(out)


def _strings(body: str):
    """Split on top-level commas; adjacent literals concatenate like in C."""
    body = _strip_comments(body)
    items = []
    for chunk in re.split(r',(?=(?:[^"\\]|\\.|"(?:[^"\\]|\\.)*")*$)', body):
        lits = _STR.findall(chunk)
        if lits:
            items.append(''.join(lits).encode().decode('unicode_escape'))
        elif chunk.strip():
            raise ValueError(f"unexpected token in string array: {chunk.strip()!r}")
    return items


def _ints(body: str):
    return [int(x) for x in re.findall(r'-?\d+', _strip_comments(body))]


def classify(help_text: str, needs_arg: int) -> str:
    """Infer a value type from the help-string prefix."""
    h = help_text.lower()
    if needs_arg == 0:
        return "flag"
    if h.startswith("(integer"):
        return "int64" if "64-bit" in h else "int"
    if h.startswith("(floating"):
        return "float"
    if h.startswith("(start,stop)"):
        return "range"
    return "str"


def parse(cmdoptions_c: str) -> list:
    src = Path(cmdoptions_c).read_text(encoding="utf-8", errors="replace")
    names = _strings(_array_body(src, "OptionArray"))
    helps = _strings(_array_body(src, "OptionHelp"))
    needs = _ints(_array_body(src, "needsArg"))
    if not (len(names) == len(helps) == len(needs)):
        raise ValueError(f"array length mismatch: names={len(names)} "
                         f"help={len(helps)} needsArg={len(needs)}")
    opts = []
    for name, h, n in zip(names, helps, needs):
        # "(Integer < 32-bit): text"  ->  text
        text = h.split(":", 1)[1].strip() if ":" in h[:22] else h.strip()
        opts.append({"name": name, "help": text, "help_raw": h,
                     "needs_arg": n, "type": classify(h, n)})
    return opts


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cmdoptions_c")
    ap.add_argument("-o", "--out", default=str(Path(__file__).with_name("yafu_options.json")))
    a = ap.parse_args()
    opts = parse(a.cmdoptions_c)
    Path(a.out).write_text(json.dumps(opts, indent=1))
    print(f"wrote {len(opts)} options to {a.out}")
