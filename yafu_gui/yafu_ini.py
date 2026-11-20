"""
yafu_ini.py -- read and write yafu.ini without disturbing its comments.

Parsing mirrors readINI() in cmdOptions.c exactly:
  * a line whose first char is '%' or ' ' is skipped
  * trailing CR/LF stripped, blank lines skipped
  * key = text before the first '=', value = text up to the NEXT '='
    (strtok semantics: "a=b=c" gives value "b"); no '=' -> flag, value None
  * repeated keys are all applied in order (tune_info relies on this)

Writing never reflows the file.  Each option edit is applied line-by-line:
  enable/change -> rewrite the active "key=..." line, else un-comment the first
                   "% key..." example line, else append under a GUI section
  disable       -> comment the active line(s) out with "% "
"""
from pathlib import Path
import re
import shutil

GUI_SECTION = "% ---- added by yafu GUI ----"
MULTI_KEYS = {"tune_info", "v"}          # keys that legitimately repeat


def _parse_active(line: str):
    """Return (key, value) as readINI would, or None if readINI skips it."""
    s = line.rstrip("\r\n")
    if not s or s[0] in "% ":
        return None
    if "=" in s:
        key, _, rest = s.partition("=")
        val = rest.split("=", 1)[0]
        return key, val
    return s, None


_COMMENTED = re.compile(r"^%\s*([A-Za-z_][A-Za-z0-9_]*)\s*(=.*)?$")


class YafuIni:
    def __init__(self, path):
        self.path = Path(path)
        self.lines = []
        if self.path.exists():
            self.lines = self.path.read_text(encoding="utf-8", errors="replace").splitlines()

    # ---------------------------------------------------------------- read
    def entries(self):
        """All active (key, value) pairs, in file order."""
        out = []
        for ln in self.lines:
            kv = _parse_active(ln)
            if kv:
                out.append(kv)
        return out

    def values(self):
        """Effective value per key (last one wins; 'v' counted; flags -> True)."""
        d = {}
        for k, v in self.entries():
            if k == "v":
                d["v"] = d.get("v", 0) + 1
            elif k in MULTI_KEYS:
                d.setdefault(k, []).append(v)
            else:
                d[k] = True if v is None else v
        return d

    # --------------------------------------------------------------- write
    def _active_idx(self, key):
        return [i for i, ln in enumerate(self.lines)
                if (kv := _parse_active(ln)) and kv[0] == key]

    def set(self, key, value):
        """value None -> flag line 'key'; otherwise 'key=value'."""
        text = key if value is None else f"{key}={value}"
        idx = self._active_idx(key)
        if idx:
            self.lines[idx[0]] = text
            for i in reversed(idx[1:]):           # collapse duplicates
                self.lines[i] = "% " + self.lines[i]
            return
        for i, ln in enumerate(self.lines):     # un-comment an example line
            m = _COMMENTED.match(ln.strip())
            if m and m.group(1) == key:
                self.lines[i] = text
                return
        if GUI_SECTION not in self.lines:
            self.lines += ["", GUI_SECTION]
        self.lines.append(text)

    def unset(self, key):
        for i in self._active_idx(key):
            self.lines[i] = "% " + self.lines[i]

    def set_verbosity(self, n):
        """readINI applies each bare 'v' line once, so level N = N lines."""
        self.unset("v")
        if n <= 0:
            return
        self.set("v", None)                      # positions the first one
        i = self._active_idx("v")[0]
        self.lines[i + 1:i + 1] = ["v"] * (n - 1)

    def save(self, backup=True):
        if backup and self.path.exists():
            shutil.copy2(self.path, self.path.with_suffix(".ini.bak"))
        self.path.write_text("\n".join(self.lines) + "\n", encoding="utf-8")
