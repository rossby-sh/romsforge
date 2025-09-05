# log_utils.py — pretty, pinpointing logs with optional ANSI colors
# API:
#   from log_utils import (
#       configure, step, capture_warnings,
#       info, note, plus, warn_line, done, ellipsis, bar
#   )
#
# Defaults:
#   - width = auto (80 if undetectable)
#   - color_mode = 'auto'  # 'on' | 'off' | 'auto'
#   - show_sections = False (no "-- section --" rulers)
#
# Behavior:
#   - In 'auto' mode, colors are enabled only when sys.stdout.isatty() is True,
#     and disabled otherwise (e.g., when redirecting to a file).
#   - You can force colors ON via configure(color_mode='on'), or OFF via
#     configure(color_mode='off'). Environment variables respected:
#       * FORCE_COLOR=1 -> treat as color_mode='on'
#       * NO_COLOR=1    -> treat as color_mode='off'

from __future__ import annotations
import os
import sys
import time
import warnings
import traceback
import contextlib
from typing import Optional
import re
import shutil

# -----------------------------------------------------------------------------
# Configuration state
TERMW = 80
SHOW_SECTIONS = False
COLOR_MODE: str = 'auto'  # 'on' | 'off' | 'auto'

# ANSI codes
RESET = "[0m"
GRAY  = "[90m"
RED   = "[91m"
BLUE  = "[94m"
YELLOW= "[93m"
GREEN = "[92m"

ANSI_RE = re.compile(r"\[[0-9;]*m")

# -----------------------------------------------------------------------------
# Helpers

def _detect_term_width() -> int:
    try:
        w = shutil.get_terminal_size().columns
        return max(40, min(160, w))
    except Exception:
        return TERMW


def _colors_enabled() -> bool:
    # Env overrides
    if os.getenv('NO_COLOR', '').strip():
        return False
    if os.getenv('FORCE_COLOR', '').strip():
        return True
    # Config overrides
    if COLOR_MODE == 'on':
        return True
    if COLOR_MODE == 'off':
        return False
    # auto
    try:
        return sys.stdout.isatty()
    except Exception:
        return False


def _c(code: str, s: str) -> str:
    if _colors_enabled():
        return f"{code}{s}{RESET}"
    return s


def strip_ansi(s: str) -> str:
    return ANSI_RE.sub('', s)

# -----------------------------------------------------------------------------
# Public configuration

def configure(width: Optional[int] = None, show_sections: Optional[bool] = None, color_mode: Optional[str] = None):
    global TERMW, SHOW_SECTIONS, COLOR_MODE
    if width is None:
        TERMW = _detect_term_width()
    else:
        TERMW = int(width)
    if show_sections is not None:
        SHOW_SECTIONS = bool(show_sections)
    if color_mode is not None:
        if color_mode not in {'on','off','auto'}:
            raise ValueError("color_mode must be 'on' | 'off' | 'auto'")
        COLOR_MODE = color_mode

# -----------------------------------------------------------------------------
# Visual helpers

def _ruler(prefix: str, title: str) -> str:
    # Build uncolored line to compute padding, then add color.
    base = f"{prefix} {title} "
    pad_len = max(0, TERMW - len(base))
    line = base + ("-" * pad_len)
    return _c(GRAY, line)


def bar(title: str):
    # Top-level banner
    base = f"== {title} "
    pad_len = max(0, TERMW - len(base))
    line = base + ("=" * pad_len)
    print(_c(GRAY, line))


def ellipsis(path: str, keep: int = 2, full_path: bool = False) -> str:
    if full_path:
        return path
    parts = path.split(os.sep)
    if len(parts) <= keep + 1:
        return path
    return os.sep + "..." + os.sep + os.sep.join(parts[-keep:])

# -----------------------------------------------------------------------------
# One-liners (status)

def info(msg: str):
    print(_c(GRAY, f"· {msg}"))


def note(msg: str):
    print(f"{_c(RED, '[NOTE]')} {msg}")


def plus(msg: str):
    print(f"{_c(GREEN, '[+]')} {msg}")


def warn_line(msg: str):
    print(f"{_c(YELLOW, '[WARN]')} {msg}")


def done(ts: str, dur: float):
    # Example: [DONE] 2025-05-01 00 |  3.454s
    print(f"{_c(BLUE, '[DONE]')} {ts:<13} | {dur:>6.3f}s")

# -----------------------------------------------------------------------------
# Step context
@contextlib.contextmanager
def step(title: str, **kv):
    """Pretty step block.
    - No section rulers unless SHOW_SECTIONS=True
    - Prints [OK]/[FAIL] with duration; re-raises exceptions (stack preserved)
    """
    meta = " ".join(f"{k}={v}" for k, v in kv.items()) if kv else ""
    head = (title if not meta else f"{title} | {meta}")
    if SHOW_SECTIONS:
        print(_ruler('--', head))
    t0 = time.time()
    try:
        yield
    except Exception as e:
        el = time.time() - t0
        print(f"{_c(RED, '[FAIL]')} {title} | dur={el:.3f}s | {type(e).__name__}: {e}")
        tb = traceback.extract_tb(e.__traceback__)
        if tb:
            src = tb[-1]
            print(f"       at {os.path.basename(src.filename)}:{src.lineno} in {src.name}")
        raise
    else:
        el = time.time() - t0
        print(f"{_c(GREEN, '[OK]')}   {title} | dur={el:.3f}s")

# -----------------------------------------------------------------------------
# Warnings capture
@contextlib.contextmanager
def capture_warnings(tag: Optional[str] = None):
    """Collect numpy/python warnings and flush as pretty [WARN] lines."""
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        yield
        for wi in w:
            fn = os.path.basename(getattr(wi, 'filename', '?'))
            ln = getattr(wi, 'lineno', '?')
            msg = str(wi.message)
            lbl = f"{fn}:{ln}"
            warn_line(f"{lbl} — {msg}" + (f" | {tag}" if tag else ''))

# Try to enable numpy warnings so capture_warnings sees them
try:
    import numpy as _np
    _np.seterr(all='warn')
except Exception:
    pass

# Initialize defaults
configure()  # sets TERMW from terminal if available, keeps other defaults

