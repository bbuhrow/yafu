#!/usr/bin/env python3
"""
yafu_gui.py -- a calculator-style desktop front end for YAFU (prototype).

    pip install PySide6
    python yafu_gui.py [--yafu path/to/yafu] [--cwd dir]

Design
------
* YAFU itself is unchanged: every job is one `yafu "<expr>" -opt val ...`
  process, exactly what a user would type.  The GUI is a command builder plus
  a live console.
* Options come from yafu_options.json, generated from cmdOptions.c by
  yafu_schema.py, so help text / arg-ness never drift from the C source.
  MODE_GROUPS below only decides *where* options appear.
* Precedence matches yafu: C defaults (initOpt) < yafu.ini < command line.
  The GUI shows the effective default in each field's placeholder and passes
  ONLY the options you change, highlighted in the panel and visible in the
  command preview.
* POSIX: yafu runs on a pseudo-terminal.  driver.c treats a non-TTY stdin as a
  batchfile and a piped stdout is block-buffered, so a plain pipe breaks both
  interactive-style runs and live progress.  Windows: QProcess pipes; see
  README for the optional 6-line driver.c patch (YAFU_GUI env var).
"""
import ast
import codecs
import json
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

from PySide6.QtCore import (QObject, QProcess, QSettings, QSocketNotifier, Qt,
                            QTimer, Signal, QRegularExpression, QProcessEnvironment)
from PySide6.QtGui import (QAction, QFontDatabase, QKeySequence,
                           QRegularExpressionValidator, QTextCursor, QColor)
from PySide6.QtWidgets import (
    QApplication, QButtonGroup, QCheckBox, QComboBox, QDialog, QDialogButtonBox,
    QFileDialog, QFormLayout, QFrame, QGridLayout, QHBoxLayout, QHeaderView, QLabel,
    QLineEdit, QListWidget, QListWidgetItem, QMainWindow, QMessageBox, QPlainTextEdit,
    QPushButton, QScrollArea, QSizePolicy, QSpinBox, QSplitter, QTabWidget,
    QTableWidget, QTableWidgetItem, QToolButton, QTreeWidget, QTreeWidgetItem,
    QVBoxLayout, QWidget)

from yafu_ini import YafuIni

APP_DIR = Path(__file__).resolve().parent
SCHEMA = {o["name"]: o for o in json.loads((APP_DIR / "yafu_options.json").read_text())}
IS_WIN = sys.platform.startswith("win")

# ---------------------------------------------------------------------------
# Defaults from initOpt() in cmdOptions.c (display only; yafu applies its own).
# Options whose C default is 0/"" meaning "auto / use tables" are omitted, so
# their field shows "auto".
C_DEFAULTS = {
    "threads": "1", "logfile": "factor.log", "session": "session.log", "obase": "10",
    "pretest_ratio": "0.3077", "xover": "95", "snfs_xover": "95", "plan": "normal",
    "stopbase": "10", "o": "nfs.dat", "job": "nfs.job", "lathreads": "1",
    "filt_bump": "5", "testsieve": "160", "psearch": "avg", "pbatch": "5000",
    "ggnfs_dir": ".\\" if IS_WIN else "./", "poly_percent_max": "8",
    "poly_testsieve": "390", "soe_block": "32768", "aprcl_p": "500", "aprcl_d": "200",
    "nprp": "1", "qssave": "siqs.dat", "siqsSSalloc": "12", "siqsMFBQ": "3.3",
    "siqsBDiv": "3.0", "siqsBT": "1000000", "inmem": "70", "keeppoly": "3",
    "avxecm_curve": "2", "B1pm1": "100000", "B1pp1": "20000", "B1ecm": "11000",
    "rhomax": "1000", "fmtmax": "1000000", "vpp1_work_file": "pp1_work.ini",
    "vpm1_work_file": "pm1_work.ini", "jsonlog": "factor.json",
    "bayes_prior": "random", "siqsMid": "auto", "siqsLarge": "auto",
    "siqsKernel": "auto", "siqsSparseKernel": "auto", "siqsLPisa": "best",
}

CHOICES = {
    "plan": ["normal", "light", "deep", "none", "noecm", "custom", "bayesian"],
    "psearch": ["avg", "min", "good", "fast", "deep", "wide"],
    "obase": ["10", "16", "8"],
    "siqsMid": ["bucket", "sparse"], "siqsLarge": ["bucket", "direct"],
    "siqsKernel": ["bitloop", "sink", "creg", "cmem"],
    "siqsSparseKernel": ["bitloop", "creg", "cmem"],
    "siqsLPisa": ["scalar", "sse41", "avx2", "avx512", "neon"],
    "siever": ["11", "12", "13", "14", "15", "16"],
    "avxecm_stg2": ["0", "1", "2"], "avxecm_curve": ["0", "1", "2"],
    "bayes_prior": ["random", "hard"],
}
DIR_OPTS = {"ggnfs_dir", "cado_dir"}
FILE_OPTS = {"ecm_path", "convert_poly_path", "nfs_params", "cuda_sieve", "job", "resume"}
SCI_OPTS = {"B1ecm"}          # applyOpt accepts 1e6-style only for B1ecm

LABELS = {"v": "Verbosity", "p": "Idle priority", "o": "NFS output file",
          "a": "Algebraic-side sieving", "r": "Rational-side sieving",
          "R": "Resume NFS job", "e": "Expression", "one": "Stop after one factor"}

# Groups shown in the options panel.  Every option in the schema that is not
# placed here (and not in HANDLED_BY_MODE) lands in "Other options" automatically.
MODE_GROUPS = {
    "General": ["threads", "v", "silent", "obase", "terse", "logfile", "session",
                "seed", "p", "no_clk_test", "vproc", "jsonlog", "jsonpretty"],
    "Factor": ["plan", "pretest_ratio", "pretest", "work", "xover", "snfs_xover",
               "noecm", "max_siqs", "max_nfs", "bayes_prior"],
    "Stop conditions": ["one", "stopk", "stop_strict", "stoplt", "stople", "stopeq",
                        "stopgt", "stopge", "stopbase", "stopprime"],
    "Output files": ["op", "of", "ou", "no_expr"],
    "SIQS": ["siqsB", "siqsTF", "siqsR", "siqsT", "siqsNB", "siqsM", "siqsLPB",
             "forceDLP", "forceTLP", "forceQLP", "inmem", "qssave", "noopt"],
    "SIQS advanced": ["siqsMFBD", "siqsMFBT", "siqsMFBQ", "siqsBDiv", "siqsBT",
                      "siqsNobat", "siqsTFSm", "siqsSPV", "siqsSSidx", "siqsSSalloc",
                      "siqsMid", "siqsLarge", "siqsKernel", "siqsSparseKernel",
                      "siqsLPisa", "noneon"],
    "NFS": ["ggnfs_dir", "gnfs", "siever", "lathreads", "o", "job", "ggnfsT",
            "a", "r", "R", "minrels", "filt_bump", "td", "nfs_batch_3lp", "keep_afb"],
    "NFS poly select": ["psearch", "pbatch", "poly_testsieve", "poly_percent_max",
                        "degree", "mindegree", "maxdegree", "keeppoly", "testsieve",
                        "skipSNFScheck", "nfs_params", "nfs_stage1_args", "keeppoly"],
    "NFS stages": ["np", "np1", "nps", "npr", "ns", "nc", "nc1", "nc2", "nc3", "ncr"],
    "NFS external tools": ["cadoMsieve", "cado_dir", "convert_poly_path",
                           "cuda_sieve", "cuda_dev"],
    "ECM": ["B1ecm", "B2ecm", "sigma", "ecm_path", "ext_ecm", "prefer_gmpecm",
            "prefer_gmpecm_stg2", "prefer_avxecm_stg2", "saveB1"],
    "AVX-ECM / GPU": ["avxecm_rns_bits", "avxecm_stg2", "avxecm_stg2mem", "avxecm_bs",
                      "avxecm_curve", "use_gpuecm", "use_gpudev", "gpucurves", "cgbn"],
    "P-1 / P+1": ["B1pm1", "B2pm1", "B1pp1", "B2pp1", "vpm1_work_file",
                  "vpp1_work_file", "resume"],
    "Rho / Fermat": ["rhomax", "fmtmax"],
    "Primality": ["nprp", "aprcl_p", "aprcl_d"],
    "Sieve of Eratosthenes": ["soe_block"],
}
# Options driven by the mode's own controls (or by the GUI itself).
# tune_info is machine-written by tune() and left alone.
HANDLED_BY_MODE = {"pscreen", "pfile", "analysis", "batchfile", "script", "e", "tune_info"}
OTHER = "Other options"

ECM_B1_PRESETS = [("t20", "11000"), ("t25", "50000"), ("t30", "250000"),
                  ("t35", "1000000"), ("t40", "3000000"), ("t45", "11000000"),
                  ("t50", "43000000"), ("t55", "110000000"), ("t60", "260000000"),
                  ("t65", "850000000")]


@dataclass
class Mode:
    key: str
    label: str
    groups: list                 # expanded by default
    extra: list = field(default_factory=list)    # shown collapsed
    hint: str = ""


ALL_ALGO = ["ECM", "AVX-ECM / GPU", "P-1 / P+1", "Rho / Fermat", "SIQS",
            "SIQS advanced", "NFS", "NFS poly select", "NFS stages",
            "NFS external tools", "Primality"]
MODES = [
    Mode("factor", "factor()", ["General", "Factor", "Stop conditions"],
         ["Output files"] + ALL_ALGO + [OTHER],
         "Number or expression, e.g. fib(300), or randb(280). Several lines = batch."),
    Mode("siqs", "SIQS", ["General", "SIQS"], ["SIQS advanced", OTHER],
         "Composite to factor with the self-initializing quadratic sieve."),
    Mode("nfs", "NFS", ["General", "NFS", "NFS poly select"],
         ["NFS stages", "NFS external tools", OTHER],
         "Composite (or SNFS form) to factor with the number field sieve."),
    Mode("ecm", "ECM", ["General", "ECM"], ["AVX-ECM / GPU", OTHER],
         "Number to run ECM curves on."),
    Mode("pm1", "P±1", ["General", "P-1 / P+1"], [OTHER], "Number for P-1 / P+1."),
    Mode("rho", "Rho / Fermat", ["General", "Rho / Fermat"], [OTHER],
         "Number for Pollard rho, Fermat, or trial division."),
    Mode("soe", "Primes (SoE)", ["General", "Sieve of Eratosthenes"], [OTHER],
         "Set the range below."),
    Mode("prime", "Primality", ["General", "Primality"], [OTHER],
         "Number to test."),
    Mode("calc", "Calculator", ["General"], ["Factor"] + ALL_ALGO + [OTHER],
         "Any yafu expression. Several lines run as a script."),
]


# ---------------------------------------------------------------------------
def quick_size(expr: str):
    """Evaluate simple integer arithmetic locally to show C-size / bits.
    Returns None when the expression uses yafu functions we don't mirror."""
    s = expr.strip().replace("^", "**")
    if not s:
        return None
    try:
        tree = ast.parse(s, mode="eval")
    except SyntaxError:
        return None

    def ev(n):
        if isinstance(n, ast.Expression):
            return ev(n.body)
        if isinstance(n, ast.Constant) and isinstance(n.value, int):
            return n.value
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.USub, ast.UAdd)):
            v = ev(n.operand)
            return -v if isinstance(n.op, ast.USub) else v
        if isinstance(n, ast.BinOp):
            a, b = ev(n.left), ev(n.right)
            if isinstance(n.op, ast.Add): return a + b
            if isinstance(n.op, ast.Sub): return a - b
            if isinstance(n.op, ast.Mult): return a * b
            if isinstance(n.op, (ast.Div, ast.FloorDiv)): return a // b
            if isinstance(n.op, ast.Mod): return a % b
            if isinstance(n.op, ast.Pow):
                if b < 0 or abs(a).bit_length() * b > 400000:
                    raise ValueError("too big to preview")
                return a ** b
        raise ValueError("unsupported")
    try:
        v = ev(tree)
    except Exception:
        return None
    return abs(v)


ANSI = re.compile(r"\x1b\[[0-9;?]*[A-Za-z]")


class TermView(QPlainTextEdit):
    """Console that honours \\r (overwrite line) and \\b like a terminal."""
    lineCompleted = Signal(str)

    def __init__(self):
        super().__init__()
        self.setReadOnly(True)
        self.setLineWrapMode(QPlainTextEdit.NoWrap)
        self.setMaximumBlockCount(200_000)
        f = QFontDatabase.systemFont(QFontDatabase.FixedFont)
        f.setPointSize(10)
        self.setFont(f)
        self.setObjectName("term")
        self._cur, self._col = "", 0

    def clear_all(self):
        self.clear()
        self._cur, self._col = "", 0

    def feed(self, text: str):
        text = ANSI.sub("", text)
        sb = self.verticalScrollBar()
        at_bottom = sb.value() >= sb.maximum() - 4
        done = []
        for piece in re.split(r"(\r\n|\n|\r|\x08)", text):
            if piece in ("\n", "\r\n"):
                done.append(self._cur)
                self._cur, self._col = "", 0
            elif piece == "\r":
                self._col = 0
            elif piece == "\x08":
                self._col = max(0, self._col - 1)
            elif piece:
                c = self._col
                self._cur = self._cur[:c] + piece + self._cur[c + len(piece):]
                self._col = c + len(piece)
        cur = self.textCursor()
        cur.movePosition(QTextCursor.End)
        cur.movePosition(QTextCursor.StartOfBlock, QTextCursor.KeepAnchor)
        cur.insertText("\n".join(done + [self._cur]))
        if at_bottom:
            sb.setValue(sb.maximum())
        for ln in done:
            self.lineCompleted.emit(ln)

    def note(self, text):
        """GUI-originated message on its own line."""
        if self._cur:
            self.feed("\n")
        self.feed(text + "\n")


# ---------------------------------------------------------------------------
class Runner(QObject):
    output = Signal(str)
    finished = Signal(int)

    def start(self, argv, cwd, env): ...
    def interrupt(self): ...
    def kill(self): ...


class PtyRunner(Runner):
    """POSIX: run yafu on a pseudo-terminal so isatty() is true and stdout is
    line-buffered -- yafu then behaves exactly as in a shell."""

    def start(self, argv, cwd, env):
        import pty, termios
        self._dec = codecs.getincrementaldecoder("utf-8")("replace")
        master, slave = pty.openpty()
        attrs = termios.tcgetattr(slave)
        attrs[1] &= ~termios.ONLCR                 # no \n -> \r\n translation
        attrs[3] &= ~(termios.ECHO | termios.ICANON)
        termios.tcsetattr(slave, termios.TCSANOW, attrs)
        try:
            self.proc = subprocess.Popen(argv, cwd=cwd, env=env, stdin=slave,
                                         stdout=slave, stderr=slave,
                                         start_new_session=True, close_fds=True)
        finally:
            os.close(slave)
        self.master = master
        self.notifier = QSocketNotifier(master, QSocketNotifier.Read, self)
        self.notifier.activated.connect(self._read)

    def _read(self):
        try:
            data = os.read(self.master, 65536)
        except OSError:
            data = b""
        if data:
            self.output.emit(self._dec.decode(data))
            return
        self.notifier.setEnabled(False)
        os.close(self.master)
        self._poll = QTimer(self, interval=50)
        self._poll.timeout.connect(self._reap)
        self._poll.start()

    def _reap(self):
        rc = self.proc.poll()
        if rc is not None:
            self._poll.stop()
            self.finished.emit(rc)

    def _signal(self, sig):
        try:
            os.killpg(self.proc.pid, sig)          # reaches child sievers too
        except (ProcessLookupError, PermissionError):
            pass

    def interrupt(self): self._signal(signal.SIGINT)
    def kill(self): self._signal(signal.SIGKILL)


class PipeRunner(Runner):
    """QProcess with pipes (Windows, or POSIX when yafu has the YAFU_GUI patch)."""

    def start(self, argv, cwd, env):
        self._dec = codecs.getincrementaldecoder("utf-8")("replace")
        self.p = QProcess(self)
        self.p.setProcessChannelMode(QProcess.MergedChannels)
        self.p.setStandardInputFile(QProcess.nullDevice())
        qenv = QProcessEnvironment()
        for k, v in env.items():
            qenv.insert(k, v)
        self.p.setProcessEnvironment(qenv)
        self.p.setWorkingDirectory(cwd)
        self.p.readyReadStandardOutput.connect(
            lambda: self.output.emit(self._dec.decode(bytes(self.p.readAllStandardOutput()))))
        self.p.finished.connect(lambda code, _st: self.finished.emit(code))
        self.p.errorOccurred.connect(self._err)
        self.p.start(argv[0], argv[1:])

    def _err(self, e):
        if e == QProcess.FailedToStart:
            self.output.emit(f"\n[failed to start: {self.p.errorString()}]\n")
            self.finished.emit(-1)

    def interrupt(self):
        if IS_WIN:      # no portable console Ctrl-C to a child; see README
            self.p.kill()
        else:
            os.kill(int(self.p.processId()), signal.SIGINT)

    def kill(self): self.p.kill()


# ---------------------------------------------------------------------------
class OptionRow:
    """One option: label + editor + reset.  `value()` is what the user set;
    `args()` is what goes on the command line (empty unless changed)."""

    def __init__(self, name, panel):
        self.name, self.panel = name, panel
        self.spec = SCHEMA.get(name, {"name": name, "help": "", "type": "str", "needs_arg": 1})
        self.type = self.spec["type"]
        self.baseline = ""          # effective default as a string, "" = auto
        self.ini_flag = False
        self.label = QLabel(LABELS.get(name, name))
        self.label.setMinimumWidth(130)
        tip = f"<b>-{name}</b><br>{self.spec.get('help_raw', '').strip()}"
        self.label.setToolTip(tip)
        self.reset_btn = QToolButton(text="↺", toolTip="Reset to default", autoRaise=True)
        self.reset_btn.clicked.connect(self.reset)
        self.editor = self._make_editor()
        self.editor.setToolTip(tip)

    # ---- editors
    def _make_editor(self):
        n, t = self.name, self.type
        if n == "v":
            w = QSpinBox(minimum=0, maximum=6)
            w.valueChanged.connect(self.changed)
            return w
        if t == "flag":
            w = QCheckBox()
            w.toggled.connect(self.changed)
            return w
        if self.spec.get("needs_arg") == 2:          # flag with optional argument
            box = QWidget()
            h = QHBoxLayout(box); h.setContentsMargins(0, 0, 0, 0)
            self.cb = QCheckBox(); self.le = QLineEdit(placeholderText="optional argument")
            if t == "range":
                self.le.setPlaceholderText("start,stop (optional)")
            h.addWidget(self.cb); h.addWidget(self.le, 1)
            self.cb.toggled.connect(self.changed); self.le.textChanged.connect(self.changed)
            return box
        if n in CHOICES:
            w = QComboBox()
            w.addItem("default", "")
            for c in CHOICES[n]:
                w.addItem(c, c)
            w.currentIndexChanged.connect(self.changed)
            return w
        w = QLineEdit()
        rx = {"int": r"\d*", "int64": r"\d*",
              "float": r"\d*\.?\d*([eE][-+]?\d*)?"}.get(t)
        if n in SCI_OPTS:
            rx = r"\d*(\.\d*)?([eE]\d*)?"
        if rx:
            w.setValidator(QRegularExpressionValidator(QRegularExpression(rx)))
        w.textChanged.connect(self.changed)
        if n in DIR_OPTS or n in FILE_OPTS:
            box = QWidget()
            h = QHBoxLayout(box); h.setContentsMargins(0, 0, 0, 0); h.setSpacing(2)
            b = QToolButton(text="…")
            b.clicked.connect(lambda: self._browse(w, n in DIR_OPTS))
            h.addWidget(w, 1); h.addWidget(b)
            self.le = w
            return box
        self.le = w
        return w

    def _browse(self, le, is_dir):
        start = self.panel.window().cwd
        p = (QFileDialog.getExistingDirectory(self.editor, self.name, start) if is_dir
             else QFileDialog.getOpenFileName(self.editor, self.name, start)[0])
        if p:
            if is_dir and not p.endswith(("/", "\\")):
                p += os.sep                     # yafu concatenates file names
            le.setText(p)

    # ---- state
    def set_baseline(self, ini_value):
        """ini_value: None (absent), True (flag present), str, or int (v count)."""
        n = self.name
        if n == "v":
            self.baseline = str(ini_value or 0)
            self.editor.blockSignals(True); self.editor.setValue(int(self.baseline))
            self.editor.blockSignals(False)
        elif self.type == "flag":
            self.ini_flag = ini_value is True
            self.editor.blockSignals(True)
            self.editor.setChecked(self.ini_flag)
            self.editor.setEnabled(not self.ini_flag)
            self.editor.blockSignals(False)
            if self.ini_flag:
                self.editor.setToolTip("Enabled in yafu.ini. Command-line flags can only "
                                       "turn options on; edit it under Settings ▸ Defaults.")
        else:
            src = "yafu.ini" if isinstance(ini_value, str) else "default"
            val = ini_value if isinstance(ini_value, str) else C_DEFAULTS.get(n, "")
            self.baseline = val
            ph = f"{val}  ({src})" if val else "auto"
            if isinstance(self.editor, QComboBox):
                self.editor.setItemText(0, f"{val or 'auto'}  ({src})")
            elif hasattr(self, "le"):
                if self.spec.get("needs_arg") != 2:
                    self.le.setPlaceholderText(ph)
        self.changed()

    def dirty(self):
        n, e = self.name, self.editor
        if n == "v":
            return str(e.value()) != self.baseline
        if self.type == "flag":
            return e.isChecked() and not self.ini_flag
        if self.spec.get("needs_arg") == 2:
            return self.cb.isChecked()
        if isinstance(e, QComboBox):
            return e.currentData() not in ("", self.baseline)
        t = self.le.text().strip()
        return bool(t) and t != self.baseline

    def args(self):
        if not self.dirty():
            return []
        n, e = self.name, self.editor
        if n == "v":
            k = e.value() - int(self.baseline)
            return ["-v"] * max(0, k)      # can only add verbosity via cmd line
        if self.type == "flag":
            return [f"-{n}"]
        if self.spec.get("needs_arg") == 2:
            a = self.le.text().strip()
            return [f"-{n}"] + ([a] if a else [])
        v = e.currentData() if isinstance(e, QComboBox) else self.le.text().strip()
        return [f"-{n}", v]

    def reset(self):
        e = self.editor
        if self.name == "v":
            e.setValue(int(self.baseline or 0))
        elif self.type == "flag":
            e.setChecked(self.ini_flag)
        elif self.spec.get("needs_arg") == 2:
            self.cb.setChecked(False); self.le.clear()
        elif isinstance(e, QComboBox):
            e.setCurrentIndex(0)
        else:
            self.le.clear()

    def changed(self, *_):
        d = self.dirty()
        f = self.label.font(); f.setBold(d); self.label.setFont(f)
        self.label.setProperty("dirty", d)
        self.label.style().unpolish(self.label); self.label.style().polish(self.label)
        self.reset_btn.setVisible(d)
        self.panel.option_changed()

    def matches(self, text):
        t = text.lower()
        return (not t) or t in self.name.lower() or t in self.spec.get("help", "").lower() \
            or t in self.label.text().lower()


class CollapsibleGroup(QWidget):
    def __init__(self, title):
        super().__init__()
        self.title = title
        v = QVBoxLayout(self); v.setContentsMargins(0, 0, 0, 0); v.setSpacing(0)
        self.head = QToolButton(checkable=True, checked=True)
        self.head.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.head.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.head.setObjectName("grouphead")
        self.head.toggled.connect(self._toggle)
        self.body = QWidget()
        self.grid = QGridLayout(self.body)
        self.grid.setContentsMargins(14, 4, 4, 8); self.grid.setVerticalSpacing(3)
        self.grid.setColumnStretch(1, 1)
        v.addWidget(self.head); v.addWidget(self.body)
        self.rows = []
        self._toggle(True)

    def add(self, row: OptionRow):
        r = len(self.rows)
        self.grid.addWidget(row.label, r, 0)
        self.grid.addWidget(row.editor, r, 1)
        self.grid.addWidget(row.reset_btn, r, 2)
        self.rows.append(row)

    def _toggle(self, on):
        self.head.setArrowType(Qt.DownArrow if on else Qt.RightArrow)
        self.body.setVisible(on)

    def set_expanded(self, on):
        self.head.setChecked(on)

    def refresh_title(self):
        n = sum(r.dirty() for r in self.rows)
        self.head.setText(self.title + (f"   •  {n} changed" if n else ""))

    def filter(self, text):
        any_vis = False
        for r in self.rows:
            vis = r.matches(text)
            for w in (r.label, r.editor):
                w.setVisible(vis)
            r.reset_btn.setVisible(vis and r.dirty())
            any_vis |= vis
        return any_vis


class OptionsPanel(QWidget):
    changed = Signal()

    def __init__(self):
        super().__init__()
        v = QVBoxLayout(self); v.setContentsMargins(0, 0, 0, 0)
        top = QHBoxLayout()
        self.search = QLineEdit(placeholderText="Filter options (name or help text)")
        self.search.setClearButtonEnabled(True)
        self.search.textChanged.connect(self.apply_filter)
        reset_all = QPushButton("Reset all")
        reset_all.clicked.connect(self.reset_all)
        top.addWidget(self.search, 1); top.addWidget(reset_all)
        v.addLayout(top)

        self.scroll = QScrollArea(widgetResizable=True)
        self.scroll.setFrameShape(QFrame.NoFrame)
        inner = QWidget(); self.vbox = QVBoxLayout(inner)
        self.vbox.setContentsMargins(0, 0, 6, 0); self.vbox.setSpacing(2)
        self.scroll.setWidget(inner)
        v.addWidget(self.scroll, 1)

        extra = QHBoxLayout()
        extra.addWidget(QLabel("Extra args"))
        self.extra = QLineEdit(placeholderText="anything else, e.g. -siqsB 9000")
        self.extra.textChanged.connect(self.option_changed)
        extra.addWidget(self.extra, 1)
        v.addLayout(extra)

        self.groups, self.rows = {}, {}
        placed = set()
        for g, names in MODE_GROUPS.items():
            grp = CollapsibleGroup(g)
            for n in dict.fromkeys(names):            # de-dup, keep order
                if n in SCHEMA and n not in placed:
                    row = OptionRow(n, self); grp.add(row)
                    self.rows[n] = row; placed.add(n)
            self.groups[g] = grp
            self.vbox.addWidget(grp)
        other = CollapsibleGroup(OTHER)
        for n in SCHEMA:
            if n not in placed and n not in HANDLED_BY_MODE:
                row = OptionRow(n, self); other.add(row); self.rows[n] = row
        self.groups[OTHER] = other
        self.vbox.addWidget(other)
        self.vbox.addStretch(1)
        self._mode = MODES[0]
        self._building = False

    def set_mode(self, mode: Mode):
        self._mode = mode
        for g, w in self.groups.items():
            shown = g in mode.groups or g in mode.extra
            w.setVisible(shown)
            if shown:
                w.set_expanded(g in mode.groups or any(r.dirty() for r in w.rows))
        self.apply_filter(self.search.text())

    def visible_groups(self):
        return [self.groups[g] for g in (self._mode.groups + self._mode.extra)
                if g in self.groups]

    def apply_filter(self, text):
        for g in self.visible_groups():
            hit = g.filter(text)
            g.setVisible(hit or not text)
            if text and hit:
                g.set_expanded(True)

    def load_baselines(self, ini_values):
        self._building = True
        for n, row in self.rows.items():
            row.set_baseline(ini_values.get(n))
        self._building = False
        self.option_changed()

    def reset_all(self):
        self._building = True
        for r in self.rows.values():
            r.reset()
        self.extra.clear()
        self._building = False
        self.option_changed()

    def option_changed(self):
        if self._building:
            return
        for g in self.groups.values():
            g.refresh_title()
        self.changed.emit()

    def args(self):
        out = []
        for g in self.visible_groups():
            for r in g.rows:
                out += r.args()
        try:
            out += shlex.split(self.extra.text(), posix=not IS_WIN)
        except ValueError:
            pass
        return out

    def hidden_dirty(self):
        vis = {id(g) for g in self.visible_groups()}
        return [r.name for g in self.groups.values() if id(g) not in vis
                for r in g.rows if r.dirty()]


# ---------------------------------------------------------------------------
class ModeArgs(QWidget):
    """Per-mode argument widgets (curves, range, test type ...)."""
    changed = Signal()

    def __init__(self, panel: OptionsPanel):
        super().__init__()
        self.panel = panel
        self.h = QHBoxLayout(self); self.h.setContentsMargins(0, 0, 0, 0)
        self.pages = {}
        e = lambda *_: self.changed.emit()

        # ECM
        w = QWidget(); l = QHBoxLayout(w); l.setContentsMargins(0, 0, 0, 0)
        self.curves = QSpinBox(minimum=1, maximum=10_000_000, value=100)
        self.curves.valueChanged.connect(e)
        self.preset = QComboBox()
        self.preset.addItem("B1 preset…", "")
        for t, b1 in ECM_B1_PRESETS:
            self.preset.addItem(f"{t}  (B1={b1})", b1)
        self.preset.activated.connect(self._apply_preset)
        l.addWidget(QLabel("Curves")); l.addWidget(self.curves)
        l.addSpacing(12); l.addWidget(self.preset); l.addStretch(1)
        self.pages["ecm"] = w

        # P±1
        w = QWidget(); l = QHBoxLayout(w); l.setContentsMargins(0, 0, 0, 0)
        self.pm_method = QComboBox(); self.pm_method.addItems(["P-1 (pm1)", "P+1 (pp1)"])
        self.pp1_bases = QSpinBox(minimum=1, maximum=1000, value=1)
        self.pm_method.currentIndexChanged.connect(e); self.pp1_bases.valueChanged.connect(e)
        l.addWidget(QLabel("Method")); l.addWidget(self.pm_method)
        l.addWidget(QLabel("P+1 bases")); l.addWidget(self.pp1_bases); l.addStretch(1)
        self.pages["pm1"] = w

        # Rho / Fermat / trial
        w = QWidget(); l = QHBoxLayout(w); l.setContentsMargins(0, 0, 0, 0)
        self.rho_method = QComboBox(); self.rho_method.addItems(["rho", "fermat", "trial"])
        self.rho_limit = QLineEdit(placeholderText="limit (fermat: iterations, trial: bound)")
        self.rho_limit.setValidator(QRegularExpressionValidator(QRegularExpression(r"\d*")))
        self.rho_method.currentIndexChanged.connect(e); self.rho_limit.textChanged.connect(e)
        l.addWidget(QLabel("Method")); l.addWidget(self.rho_method); l.addWidget(self.rho_limit, 1)
        self.pages["rho"] = w

        # SoE
        w = QWidget(); l = QHBoxLayout(w); l.setContentsMargins(0, 0, 0, 0)
        self.lo = QLineEdit("0"); self.hi = QLineEdit("1000000000")
        for x in (self.lo, self.hi):
            x.textChanged.connect(e)
        self.soe_out = QComboBox()
        self.soe_out.addItems(["count only", "list to screen", "write primes.dat"])
        self.soe_twins = QCheckBox("twins")
        self.soe_out.currentIndexChanged.connect(e); self.soe_twins.toggled.connect(e)
        l.addWidget(QLabel("Low")); l.addWidget(self.lo, 1)
        l.addWidget(QLabel("High")); l.addWidget(self.hi, 1)
        l.addWidget(self.soe_out); l.addWidget(self.soe_twins)
        self.pages["soe"] = w

        # Primality
        w = QWidget(); l = QHBoxLayout(w); l.setContentsMargins(0, 0, 0, 0)
        self.ptest = QComboBox()
        self.ptest.addItem("isprime — PRP test", "isprime")
        self.ptest.addItem("bpsw — Baillie-PSW", "bpsw")
        self.ptest.addItem("aprcl — APR-CL proof", "aprcl")
        self.ptest.currentIndexChanged.connect(e)
        l.addWidget(QLabel("Test")); l.addWidget(self.ptest); l.addStretch(1)
        self.pages["prime"] = w

        for p in self.pages.values():
            self.h.addWidget(p); p.hide()

    def _apply_preset(self, i):
        b1 = self.preset.itemData(i)
        if b1:
            self.panel.rows["B1ecm"].le.setText(b1)
        self.preset.setCurrentIndex(0)

    def set_mode(self, key):
        for k, p in self.pages.items():
            p.setVisible(k == key)
        self.setVisible(key in self.pages)

    def wrap(self, key, x):
        """Function call for one input x (x may be '@' for batches)."""
        if key in ("factor", "siqs", "nfs"):
            return f"{key}({x})"
        if key == "ecm":
            return f"ecm({x},{self.curves.value()})"
        if key == "pm1":
            return f"pm1({x})" if self.pm_method.currentIndex() == 0 \
                else f"pp1({x},{self.pp1_bases.value()})"
        if key == "rho":
            m, lim = self.rho_method.currentText(), self.rho_limit.text().strip()
            if m == "rho":
                return f"rho({x})"
            return f"{m}({x},{lim or ('1000000' if m == 'fermat' else '10000')})"
        if key == "prime":
            return f"{self.ptest.currentData()}({x})"
        return x

    def soe_job(self):
        lo, hi = self.lo.text().strip() or "0", self.hi.text().strip() or "0"
        mode = self.soe_out.currentIndex()
        expr = f"primes({lo},{hi})" if mode == 0 else f"primes({lo},{hi},0)"
        args = [[], ["-pscreen"], ["-pfile"]][mode]
        if self.soe_twins.isChecked():
            args += ["-analysis", "2"]
        return expr, args


# ---------------------------------------------------------------------------
class ResultsView(QTreeWidget):
    """Parses yafu's console output into a per-input results tree."""
    RX_FAC = re.compile(r"^(P|PRP|C|U)(\d+)([xo]?) = (\S+)$")
    RX_ANS = re.compile(r"^(\w+) = (-?[0-9a-fA-Fx]+)$")
    RX_BATCH = re.compile(r"^=== Starting work on batchfile expression #(\d+) ===")

    def __init__(self):
        super().__init__()
        self.setHeaderLabels(["Input / factor", "Type", "Digits"])
        self.setRootIsDecorated(True)
        self.setUniformRowHeights(True)
        self.header().setStretchLastSection(False)
        self.header().setSectionResizeMode(0, QHeaderView.Stretch)
        self.header().setSectionResizeMode(1, QHeaderView.ResizeToContents)
        self.header().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.setSelectionMode(QTreeWidget.ExtendedSelection)
        self.cur = None
        self._expect = None

    def new_job(self, title):
        self.cur = QTreeWidgetItem([title, "", ""])
        f = self.cur.font(0); f.setBold(True); self.cur.setFont(0, f)
        self.addTopLevelItem(self.cur)
        self.cur.setExpanded(True)
        self.scrollToItem(self.cur)

    def line(self, s):
        s = s.strip()
        if self._expect == "batch":
            self._expect = None
            self.new_job(s)
            return
        if self._expect == "fz":
            self._expect = None
            self._add(s, "factorization", "")
            return
        if self.RX_BATCH.match(s):
            self._expect = "batch"
            return
        if s == "***factorization:***":
            self._expect = "fz"
            return
        m = self.RX_FAC.match(s)
        if m:
            self._add(m.group(4), {"P": "prime", "PRP": "prp", "C": "composite",
                                   "U": "unknown"}[m.group(1)], m.group(2))
            return
        m = self.RX_ANS.match(s)
        if m and m.group(1) == "ans":
            self._add(m.group(2), "ans", str(len(m.group(2).lstrip("-"))))

    def _add(self, value, kind, digits):
        if self.cur is None:
            self.new_job("(output)")
        it = QTreeWidgetItem([value, kind, digits])
        it.setToolTip(0, value)
        if kind == "prime" or kind == "prp":
            it.setForeground(1, QColor("#2e7d32"))
        elif kind == "composite":
            it.setForeground(1, QColor("#c62828"))
        self.cur.addChild(it)
        self.scrollToItem(it)

    def copy_selected(self):
        items = self.selectedItems() or ([self.currentItem()] if self.currentItem() else [])
        QApplication.clipboard().setText("\n".join(i.text(0) for i in items))


# ---------------------------------------------------------------------------
class SettingsDialog(QDialog):
    def __init__(self, win):
        super().__init__(win)
        self.setWindowTitle("Preferences")
        f = QFormLayout(self)
        self.exe = QLineEdit(win.exe)
        b1 = QToolButton(text="…")
        b1.clicked.connect(lambda: self._pick(self.exe, False))
        self.cwd = QLineEdit(win.cwd)
        b2 = QToolButton(text="…")
        b2.clicked.connect(lambda: self._pick(self.cwd, True))
        for le, b, lab in ((self.exe, b1, "yafu executable"), (self.cwd, b2, "Working directory")):
            h = QHBoxLayout(); h.addWidget(le, 1); h.addWidget(b); f.addRow(lab, h)
        self.backend = QComboBox()
        self.backend.addItems(["auto", "pty", "pipes"])
        self.backend.setCurrentText(win.backend)
        self.backend.setToolTip("auto = pseudo-terminal on Linux/macOS, pipes on Windows.\n"
                                "'pipes' on POSIX needs yafu built with the YAFU_GUI patch.")
        f.addRow("Process I/O", self.backend)
        note = QLabel("The working directory is where yafu reads yafu.ini and writes "
                      "siqs.dat, nfs.dat, factor.log, etc.")
        note.setWordWrap(True); note.setObjectName("muted")
        f.addRow(note)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        f.addRow(bb)
        self.resize(560, 0)

    def _pick(self, le, is_dir):
        p = (QFileDialog.getExistingDirectory(self, "Working directory", le.text()) if is_dir
             else QFileDialog.getOpenFileName(self, "yafu executable", le.text())[0])
        if p:
            le.setText(p)


class DefaultsDialog(QDialog):
    """Edit yafu.ini as a table: every option from the schema, enabled or not."""

    def __init__(self, win):
        super().__init__(win)
        self.win = win
        self.setWindowTitle("Defaults — yafu.ini")
        self.ini = YafuIni(Path(win.cwd) / "yafu.ini")
        vals = self.ini.values()
        v = QVBoxLayout(self)
        path_lbl = QLabel(f"{self.ini.path}" + ("" if self.ini.path.exists() else "  (will be created)"))
        path_lbl.setObjectName("muted")
        v.addWidget(path_lbl)
        self.filter = QLineEdit(placeholderText="Filter")
        self.filter.setClearButtonEnabled(True)
        self.filter.textChanged.connect(self._filter)
        v.addWidget(self.filter)
        self.t = QTableWidget(0, 6)
        self.t.setHorizontalHeaderLabels(["On", "Option", "Value", "Built-in default", "Description", "Section"])
        self.t.verticalHeader().setDefaultSectionSize(22)
        self.t.verticalHeader().hide()
        self.t.horizontalHeader().setSectionResizeMode(4, QHeaderView.Stretch)
        self.t.setColumnWidth(5, 150)
        section = {n: g for g, ns in MODE_GROUPS.items() for n in ns}
        self.t.setColumnWidth(0, 36); self.t.setColumnWidth(1, 150)
        self.t.setColumnWidth(2, 220); self.t.setColumnWidth(3, 130)
        order = [n for g in MODE_GROUPS.values() for n in g if n in SCHEMA]
        order = list(dict.fromkeys(order + [n for n in SCHEMA if n not in HANDLED_BY_MODE]))
        self._orig = {}
        for n in order:
            if n == "tune_info":
                continue            # machine-written by tune(); left untouched
            s = SCHEMA[n]
            r = self.t.rowCount(); self.t.insertRow(r)
            on = QTableWidgetItem()
            cur = vals.get(n)
            on.setCheckState(Qt.Checked if cur is not None else Qt.Unchecked)
            name = QTableWidgetItem(n); name.setFlags(name.flags() & ~Qt.ItemIsEditable)
            name.setToolTip(s["help_raw"].strip())
            if n == "v":
                val = str(cur or "")
            else:
                val = "" if cur in (None, True) else str(cur)
            val_it = QTableWidgetItem(val)
            if s["type"] == "flag" and n != "v":
                val_it.setFlags(val_it.flags() & ~Qt.ItemIsEditable)
                val_it.setText("(flag)")
            dflt = QTableWidgetItem(C_DEFAULTS.get(n, "")); dflt.setFlags(dflt.flags() & ~Qt.ItemIsEditable)
            hlp = QTableWidgetItem(s["help"]); hlp.setFlags(hlp.flags() & ~Qt.ItemIsEditable)
            hlp.setToolTip(s["help"])
            sec = QTableWidgetItem(section.get(n, OTHER)); sec.setFlags(sec.flags() & ~Qt.ItemIsEditable)
            for c, it in enumerate((on, name, val_it, dflt, hlp, sec)):
                self.t.setItem(r, c, it)
            self._orig[n] = (on.checkState() == Qt.Checked, val)
        self.t.itemChanged.connect(self._edited)
        v.addWidget(self.t, 1)
        bb = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.save); bb.rejected.connect(self.reject)
        v.addWidget(bb)
        self.resize(1000, 640)

    def _edited(self, it):
        if it.column() == 2 and it.text().strip():
            self.t.item(it.row(), 0).setCheckState(Qt.Checked)

    def _filter(self, text):
        t = text.lower()
        for r in range(self.t.rowCount()):
            hay = " ".join(self.t.item(r, c).text() for c in (1, 4, 5)).lower()
            self.t.setRowHidden(r, bool(t) and t not in hay)

    def save(self):
        changed = 0
        for r in range(self.t.rowCount()):
            n = self.t.item(r, 1).text()
            on = self.t.item(r, 0).checkState() == Qt.Checked
            flag = SCHEMA[n]["type"] == "flag" and n != "v"
            val = "" if flag else self.t.item(r, 2).text().strip()
            if (on, val) == self._orig[n]:
                continue
            changed += 1
            if n == "v":
                self.ini.set_verbosity(int(val or 1) if on else 0)
            elif not on:
                self.ini.unset(n)
            elif flag:
                self.ini.set(n, None)
            elif not val and SCHEMA[n]["needs_arg"] == 1:
                QMessageBox.warning(self, "Missing value", f"{n} needs a value.")
                return
            else:
                self.ini.set(n, val or None)
        if changed:
            self.ini.save(backup=True)
        self.accept()


# ---------------------------------------------------------------------------
class MainWindow(QMainWindow):
    def __init__(self, exe=None, cwd=None):
        super().__init__()
        self.settings = QSettings("yafu", "yafu-gui")
        self.exe = exe or self.settings.value("exe", "") or self._find_yafu()
        self.cwd = cwd or self.settings.value("cwd", "") or \
            (str(Path(self.exe).parent) if self.exe else os.getcwd())
        self.backend = self.settings.value("backend", "auto")
        self.runner = None
        self.t0 = None
        self.history = []
        self.hist_pos = None
        self.setWindowTitle("YAFU")

        # ---- left: input, mode args, preview, run bar, output/results
        self.mode_bar = QWidget()
        mb = QHBoxLayout(self.mode_bar); mb.setContentsMargins(0, 0, 0, 0); mb.setSpacing(0)
        self.mode_group = QButtonGroup(self)
        for i, m in enumerate(MODES):
            b = QToolButton(text=m.label, checkable=True)
            b.setObjectName("modebtn")
            b.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            self.mode_group.addButton(b, i); mb.addWidget(b)
        self.mode_group.idClicked.connect(self.set_mode)
        self.mode_group.button(0).setChecked(True)

        self.input = QPlainTextEdit()
        self.input.setObjectName("input")
        self.input.setFixedHeight(72)
        f = QFontDatabase.systemFont(QFontDatabase.FixedFont); f.setPointSize(12)
        self.input.setFont(f)
        self.input.textChanged.connect(self.refresh)
        self.input.installEventFilter(self)
        self.size_lbl = QLabel(); self.size_lbl.setObjectName("muted")

        self.options = OptionsPanel()
        self.options.changed.connect(self.refresh)
        self.mode_args = ModeArgs(self.options)
        self.mode_args.changed.connect(self.refresh)

        self.preview = QLineEdit(readOnly=True)
        self.preview.setObjectName("preview")
        copy_cmd = QToolButton(text="Copy")
        copy_cmd.clicked.connect(lambda: QApplication.clipboard().setText(self.preview.text()))

        self.run_btn = QPushButton("Run")
        self.run_btn.setObjectName("run")
        self.run_btn.setShortcut(QKeySequence("Ctrl+Return"))
        self.run_btn.setToolTip("Run (Ctrl+Enter)")
        self.run_btn.clicked.connect(self.run)
        self.stop_btn = QPushButton("Stop")
        self.stop_btn.setToolTip("First click sends Ctrl-C (yafu saves state where it can); "
                                 "second click kills the process.")
        self.stop_btn.clicked.connect(self.stop)
        self.stop_btn.setEnabled(False)
        clear_btn = QPushButton("Clear output")
        self.status = QLabel("Idle")

        self.term = TermView()
        self.results = ResultsView()
        self.term.lineCompleted.connect(self.results.line)
        clear_btn.clicked.connect(lambda: (self.term.clear_all(), self.results.clear()))

        res_box = QWidget(); rv = QVBoxLayout(res_box); rv.setContentsMargins(0, 0, 0, 0)
        rh = QHBoxLayout(); rh.addWidget(QLabel("Results")); rh.addStretch(1)
        cp = QToolButton(text="Copy"); cp.clicked.connect(self.results.copy_selected)
        cl = QToolButton(text="Clear"); cl.clicked.connect(self.results.clear)
        rh.addWidget(cp); rh.addWidget(cl)
        rv.addLayout(rh); rv.addWidget(self.results)

        out_split = QSplitter(Qt.Horizontal)
        out_split.addWidget(self.term); out_split.addWidget(res_box)
        out_split.setStretchFactor(0, 3); out_split.setStretchFactor(1, 2)
        out_split.setSizes([620, 420])
        self.out_split = out_split

        left = QWidget(); lv = QVBoxLayout(left); lv.setContentsMargins(8, 8, 4, 8)
        lv.addWidget(self.mode_bar)
        lv.addSpacing(4)
        lv.addWidget(self.input)
        lv.addWidget(self.size_lbl)
        lv.addWidget(self.mode_args)
        ph = QHBoxLayout(); ph.addWidget(self.preview, 1); ph.addWidget(copy_cmd)
        lv.addLayout(ph)
        bar = QHBoxLayout()
        bar.addWidget(self.run_btn); bar.addWidget(self.stop_btn); bar.addWidget(clear_btn)
        bar.addSpacing(12); bar.addWidget(self.status, 1)
        lv.addLayout(bar)
        lv.addWidget(out_split, 1)

        # ---- right: options / history
        self.hist_list = QListWidget()
        self.hist_list.itemDoubleClicked.connect(self.restore_history)
        right = QTabWidget()
        right.addTab(self.options, "Options")
        right.addTab(self.hist_list, "History")
        right.setMinimumWidth(380)

        split = QSplitter(Qt.Horizontal)
        split.addWidget(left); split.addWidget(right)
        split.setStretchFactor(0, 3); split.setStretchFactor(1, 2)
        self.split = split
        self.setCentralWidget(split)

        # ---- menu
        m = self.menuBar().addMenu("&File")
        a = QAction("Preferences…", self); a.triggered.connect(self.prefs); m.addAction(a)
        a = QAction("Open working directory", self)
        a.triggered.connect(lambda: self._open_path(self.cwd)); m.addAction(a)
        m.addSeparator()
        a = QAction("Quit", self, shortcut=QKeySequence.Quit); a.triggered.connect(self.close)
        m.addAction(a)
        m = self.menuBar().addMenu("&Settings")
        a = QAction("Defaults (yafu.ini)…", self); a.triggered.connect(self.defaults)
        m.addAction(a)
        a = QAction("Reload yafu.ini", self); a.triggered.connect(self.load_ini); m.addAction(a)
        m = self.menuBar().addMenu("&Help")
        a = QAction("yafu option list (-h)", self)
        a.triggered.connect(lambda: self.launch([self.exe, "-h"], "yafu -h"))
        m.addAction(a)
        a = QAction("Open docfile.txt", self)
        a.triggered.connect(lambda: self._open_path(str(Path(self.cwd) / "docfile.txt")))
        m.addAction(a)

        self.elapsed = QLabel()
        self.where = QLabel(); self.where.setObjectName("muted")
        self.statusBar().addWidget(self.where, 1)
        self.statusBar().addPermanentWidget(self.elapsed)
        self.tick = QTimer(self, interval=500); self.tick.timeout.connect(self._tick)

        self.load_ini()
        g = self.settings.value("geometry")
        if g is not None:
            self.restoreGeometry(g)
        else:
            self.resize(1400, 860)
        self.history = json.loads(self.settings.value("history", "[]"))
        for h in self.history:
            self._add_hist_item(h)
        mi = int(self.settings.value("mode", 0))
        self.mode_group.button(mi).setChecked(True)
        self.set_mode(mi)

    # ---------------------------------------------------------------- setup
    @staticmethod
    def _find_yafu():
        here = APP_DIR
        names = ["yafu.exe", "yafu-x64.exe"] if IS_WIN else ["yafu"]
        for d in (here, here.parent, Path.cwd()):
            for n in names:
                if (d / n).is_file():
                    return str(d / n)
        for n in names:
            if p := shutil.which(n):
                return p
        return ""

    def load_ini(self):
        ini = YafuIni(Path(self.cwd) / "yafu.ini")
        self.ini_values = ini.values()
        self.options.load_baselines(self.ini_values)
        self.where.setText(f"{self.exe or '(yafu not found — File ▸ Preferences)'}   ·   "
                           f"cwd {self.cwd}   ·   "
                           f"{'yafu.ini loaded' if ini.path.exists() else 'no yafu.ini'}")
        self.refresh()

    def prefs(self):
        d = SettingsDialog(self)
        if d.exec():
            self.exe, self.cwd = d.exe.text().strip(), d.cwd.text().strip()
            self.backend = d.backend.currentText()
            for k in ("exe", "cwd", "backend"):
                self.settings.setValue(k, getattr(self, k))
            self.load_ini()

    def defaults(self):
        if DefaultsDialog(self).exec():
            self.load_ini()
            self.term.note("[yafu.ini saved; previous version kept as yafu.ini.bak]")

    def _open_path(self, p):
        from PySide6.QtGui import QDesktopServices
        from PySide6.QtCore import QUrl
        QDesktopServices.openUrl(QUrl.fromLocalFile(p))

    # ---------------------------------------------------------------- modes
    @property
    def mode(self):
        return MODES[self.mode_group.checkedId()]

    def set_mode(self, i):
        m = MODES[i]
        self.settings.setValue("mode", i)
        self.options.set_mode(m)
        self.mode_args.set_mode(m.key)
        self.input.setPlaceholderText(m.hint)
        self.input.setEnabled(m.key != "soe")
        self.refresh()

    def eventFilter(self, obj, ev):
        # Ctrl+Up / Ctrl+Down: recall previous inputs, like the yafu prompt
        if obj is self.input and ev.type() == ev.Type.KeyPress \
                and ev.modifiers() & Qt.ControlModifier and ev.key() in (Qt.Key_Up, Qt.Key_Down):
            inputs = [h["input"] for h in self.history if h["input"]]
            if inputs:
                if self.hist_pos is None:
                    self.hist_pos = len(inputs)
                self.hist_pos += -1 if ev.key() == Qt.Key_Up else 1
                self.hist_pos = max(0, min(len(inputs) - 1, self.hist_pos))
                self.input.setPlainText(inputs[self.hist_pos])
            return True
        return super().eventFilter(obj, ev)

    # ---------------------------------------------------------------- build
    def build_job(self):
        """Return (argv, files_to_write, title).  files: {name: text} in cwd."""
        m = self.mode
        lines = [l.strip() for l in self.input.toPlainText().splitlines() if l.strip()]
        opts = self.options.args()
        files = {}
        if m.key == "soe":
            expr, extra = self.mode_args.soe_job()
            opts = extra + opts
            title = expr
        elif not lines:
            return None
        elif len(lines) > 1 and m.key == "calc":
            files[".yafu_gui_script.txt"] = "\n".join(lines) + "\n"
            expr, title = None, f"script ({len(lines)} lines)"
            opts = ["-script", ".yafu_gui_script.txt"] + opts
        elif len(lines) > 1:
            files[".yafu_gui_batch.txt"] = "\n".join(lines) + "\n"
            expr = self.mode_args.wrap(m.key, "@")
            title = f"batch: {expr} × {len(lines)}"
            opts = ["-batchfile", ".yafu_gui_batch.txt"] + opts
        else:
            expr = self.mode_args.wrap(m.key, lines[0])
            title = expr
        argv = [self.exe or "yafu"] + ([expr] if expr else []) + opts
        return argv, files, title

    def refresh(self):
        text = self.input.toPlainText()
        lines = [l for l in text.splitlines() if l.strip()]
        if self.mode.key == "soe":
            self.size_lbl.setText("")
        elif len(lines) > 1:
            self.size_lbl.setText(f"{len(lines)} inputs — runs as a "
                                  f"{'script' if self.mode.key == 'calc' else 'batchfile'}")
        elif lines:
            v = quick_size(lines[0])
            self.size_lbl.setText(f"C{len(str(v))}  ·  {v.bit_length()} bits" if v
                                  else "expression evaluated by yafu")
        else:
            self.size_lbl.setText("")
        job = self.build_job()
        if job:
            argv = job[0]
            q = subprocess.list2cmdline(argv) if IS_WIN else shlex.join(argv)
            hid = self.options.hidden_dirty()
            self.preview.setText(q)
            self.preview.setToolTip(q + (f"\n\nNot applied in this mode: {', '.join(hid)}"
                                         if hid else ""))
        else:
            self.preview.setText("")
        self.run_btn.setEnabled(bool(job) and self.runner is None)

    # ---------------------------------------------------------------- run
    def run(self):
        job = self.build_job()
        if not job or self.runner:
            return
        argv, files, title = job
        if not self.exe or not Path(self.exe).exists():
            QMessageBox.warning(self, "yafu not found",
                                "Set the yafu executable under File ▸ Preferences.")
            return
        for name, text in files.items():
            (Path(self.cwd) / name).write_text(text)
        if not title.startswith("batch"):       # batch items announce themselves
            self.results.new_job(title)
        self.launch(argv, title)
        self._record(title, argv)

    def launch(self, argv, title):
        if self.runner:
            return
        env = dict(os.environ, YAFU_GUI="1")
        use_pty = (self.backend == "pty") or (self.backend == "auto" and not IS_WIN)
        self.runner = PtyRunner(self) if use_pty else PipeRunner(self)
        self.runner.output.connect(self.term.feed)
        self.runner.finished.connect(self._done)
        self.term.note(f"\n▶ {title}")
        self.term.note("$ " + (subprocess.list2cmdline(argv) if IS_WIN else shlex.join(argv)))
        self._interrupts = 0
        try:
            self.runner.start(argv, self.cwd, env)
        except OSError as e:
            self.term.note(f"[could not start yafu: {e}]")
            self.runner = None
            return
        self.t0 = time.monotonic()
        self.tick.start()
        self._set_running(True, "Running")

    def stop(self):
        if not self.runner:
            return
        self._interrupts += 1
        if self._interrupts == 1:
            self.runner.interrupt()
            self.status.setText("Interrupt sent — waiting for yafu to save state…")
            self.stop_btn.setText("Kill")
        else:
            self.runner.kill()

    def _done(self, code):
        dt = time.monotonic() - self.t0 if self.t0 else 0
        self.tick.stop()
        self.runner.deleteLater()
        self.runner = None
        state = "Finished" if code == 0 else ("Interrupted" if self._interrupts else f"Exited ({code})")
        self.term.note(f"■ {state} in {self._fmt(dt)}")
        self._set_running(False, f"{state} in {self._fmt(dt)}")
        if self.history:
            self.history[-1]["status"] = state
            self.history[-1]["secs"] = round(dt, 1)
            self._save_hist()
            self.hist_list.item(self.hist_list.count() - 1).setText(self._hist_text(self.history[-1]))
        self.refresh()

    def _set_running(self, on, text):
        self.status.setText(text)
        self.stop_btn.setEnabled(on)
        self.stop_btn.setText("Stop")
        self.run_btn.setEnabled(not on)
        self.mode_bar.setEnabled(not on)

    def _tick(self):
        if self.t0:
            self.elapsed.setText(self._fmt(time.monotonic() - self.t0))

    @staticmethod
    def _fmt(s):
        s = int(s) if s >= 10 else s
        if isinstance(s, float):
            return f"{s:.2f} s"
        h, r = divmod(s, 3600); mi, se = divmod(r, 60)
        return f"{h}:{mi:02d}:{se:02d}" if h else f"{mi}:{se:02d}"

    # ---------------------------------------------------------------- history
    def _record(self, title, argv):
        h = {"when": time.strftime("%m-%d %H:%M"), "mode": self.mode_group.checkedId(),
             "input": self.input.toPlainText(), "title": title, "argv": argv,
             "status": "running", "secs": None}
        self.history.append(h)
        self.history = self.history[-200:]
        self.hist_pos = None
        self._add_hist_item(h)
        self._save_hist()

    def _hist_text(self, h):
        t = f"{h['secs']} s" if h.get("secs") is not None else ""
        return f"{h['when']}   {h['title'][:70]}   [{h['status']}{', ' + t if t else ''}]"

    def _add_hist_item(self, h):
        it = QListWidgetItem(self._hist_text(h))
        it.setToolTip(shlex.join(h["argv"]))
        self.hist_list.addItem(it)

    def _save_hist(self):
        self.settings.setValue("history", json.dumps(self.history))

    def restore_history(self, item):
        h = self.history[self.hist_list.row(item)]
        self.mode_group.button(h["mode"]).setChecked(True)
        self.set_mode(h["mode"])
        self.input.setPlainText(h["input"])

    def closeEvent(self, e):
        if self.runner:
            if QMessageBox.question(self, "yafu is running",
                                    "A job is still running. Kill it and quit?") != QMessageBox.Yes:
                e.ignore(); return
            self.runner.kill()
        self.settings.setValue("geometry", self.saveGeometry())
        super().closeEvent(e)


STYLE = """
QToolButton#modebtn { padding: 7px 10px; border: 1px solid palette(mid); border-left: none;
                      background: palette(button); }
QToolButton#modebtn:first { border-left: 1px solid palette(mid); }
QToolButton#modebtn:checked { background: #2f6fde; color: white; border-color: #2f6fde; }
QPlainTextEdit#input { border: 2px solid #2f6fde; border-radius: 6px; padding: 4px; }
QLineEdit#preview { font-family: monospace; color: palette(dark); background: palette(alternate-base); }
QPlainTextEdit#term { background: #1e1f22; color: #d6d6d6; border: none; }
QPushButton#run { background: #2f6fde; color: white; padding: 6px 22px; border-radius: 4px; font-weight: bold; }
QPushButton#run:disabled { background: #9db4e0; }
QToolButton#grouphead { text-align: left; padding: 5px 2px; font-weight: bold; border: none;
                        border-bottom: 1px solid palette(mid); }
QLabel[dirty="true"] { color: #2f6fde; }
QLabel#muted { color: palette(dark); }
"""


def main():
    import argparse
    ap = argparse.ArgumentParser(description="YAFU GUI prototype")
    ap.add_argument("--yafu", help="path to the yafu executable")
    ap.add_argument("--cwd", help="working directory (where yafu.ini lives)")
    a = ap.parse_args()
    app = QApplication(sys.argv[:1])
    app.setStyle("Fusion")
    app.setStyleSheet(STYLE)
    w = MainWindow(a.yafu, a.cwd)
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
