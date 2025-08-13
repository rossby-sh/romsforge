# log_utils.py — pretty, pinpointing logs (no section rulers by default)
# Usage:
#   from log_utils import step, capture_warnings, info, note, plus, warn_line, done, ellipsis, configure
#   configure(width=80, show_sections=False)  # no "-- section --" rulers

from __future__ import annotations
import os
import time
import warnings
import traceback
import contextlib
from typing import Optional

TERMW = 80
SHOW_SECTIONS = False  # when True, prints "-- title -----" rulers

RESET = "\033[0m"
GRAY  = "\033[90m"
RED   = "\033[91m"
BLUE  = "\033[94m"
YELLOW = "\033[93m"
GREEN = "\033[92m"
# --- configuration ------------------------------------------------------------
def configure(width: Optional[int] = None, show_sections: Optional[bool] = None):
    global TERMW, SHOW_SECTIONS
    if width is not None:
        TERMW = int(width)
    if show_sections is not None:
        SHOW_SECTIONS = bool(show_sections)

# --- visual helpers -----------------------------------------------------------
def _ruler(prefix: str, title: str) -> str:
    t = f"{GRAY}{prefix} {title} {RESET}"
    fill = f"{GRAY}{'-' * max(0, TERMW - len(prefix) - len(title) - 2)}{RESET}"
    return t + fill

def bar(title: str):
    # Top-level banner (rarely used; keep if you like headers at the very top)
    t = f"{GRAY}== {title} {RESET}"
    fill = f"{GRAY}{'=' * max(0, TERMW - len(title) - 4)}{RESET}"
    print(t + fill)

def ellipsis(path: str, keep: int = 2) -> str:
    parts = path.split(os.sep)
    if len(parts) <= keep + 1:
        return path
    return os.sep + "..." + os.sep + os.sep.join(parts[-keep:])

# --- one-liners ---------------------------------------------------------------

def info(msg: str):
    print(f"{GRAY}· {msg}{RESET}")

def note(msg: str):
    print(f"{RED}[NOTE]{RESET} {msg}")

def plus(msg: str):
    print(f"{GREEN}[+]{RESET} {msg}")

def warn_line(msg: str):
    print(f"{YELLOW}[WARN]{RESET} {msg}")

def done(ts: str, dur: float):
    print(f"{BLUE}[DONE]{RESET} {ts:<13} | {dur:>6.3f}s")


# --- step context -------------------------------------------------------------
@contextlib.contextmanager
def step(title: str, **kv):
    """Pretty step block. By default, NO section rulers are printed.
    Prints [OK]/[FAIL] with duration; re-raises exceptions (stack preserved)."""
    meta = " ".join(f"{k}={v}" for k, v in kv.items()) if kv else ""
    head = (title if not meta else f"{title} | {meta}")
    if SHOW_SECTIONS:
        print(_ruler("--", head))
    t0 = time.time()
    try:
        yield
    except Exception as e:
        el = time.time() - t0
        print(f"[FAIL] {title} | dur={el:.3f}s | {type(e).__name__}: {e}")
        tb = traceback.extract_tb(e.__traceback__)
        if tb:
            src = tb[-1]
            print(f"       at {os.path.basename(src.filename)}:{src.lineno} in {src.name}")
        raise
    else:
        el = time.time() - t0
        print(f"[OK]   {title} | dur={el:.3f}s")

@contextlib.contextmanager
def capture_warnings(tag: Optional[str] = None):
    """Collect numpy/python warnings and flush as pretty [WARN] lines."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        yield
        for wi in w:
            fn = os.path.basename(getattr(wi, "filename", "?"))
            ln = getattr(wi, "lineno", "?")
            msg = str(wi.message)
            lbl = f"{fn}:{ln}"
            warn_line(f"{lbl} — {msg}" + (f" | {tag}" if tag else ""))

# Make numpy raise warnings (instead of being silent) so capture_warnings sees them
try:
    import numpy as _np
    _np.seterr(all="warn")
except Exception:
    pass

