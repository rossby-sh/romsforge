#!/usr/bin/env python3
# make_concat_animation.py
import os
import glob
import re
from pathlib import Path
from PIL import Image
import numpy as np
import imageio.v2 as imageio

# -------- 유틸 --------
def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r'(\d+)', s)]

def load_images_from_dir(directory):
    files = sorted(glob.glob(os.path.join(directory, "*.png")), key=natural_key)
    if not files:
        raise FileNotFoundError(f"디렉토리에 PNG 없음: {directory}")
    return files

def resize_to_same_height(imgs, interp=Image.BICUBIC):
    h = max(im.height for im in imgs)
    out = []
    for im in imgs:
        if im.height != h:
            w = int(round(im.width * (h / im.height)))
            im = im.resize((w, h), interp)
        out.append(im)
    return out, h

def concat_h(imgs):
    imgs, h = resize_to_same_height(imgs)
    w = sum(im.width for im in imgs)
    canvas = Image.new("RGB", (w, h), (0, 0, 0))
    x = 0
    for im in imgs:
        canvas.paste(im.convert("RGB"), (x, 0))
        x += im.width
    return canvas

def ensure_even_dims(arr):
    h, w = arr.shape[:2]
    pad_h = h % 2
    pad_w = w % 2
    if pad_h or pad_w:
        new = np.zeros((h + pad_h, w + pad_w, 3), dtype=arr.dtype)
        new[:h, :w] = arr
        return new
    return arr

# -------- 메인 --------
def make_animation(left_dir, mid_dir, right_dir, out_dir, fps=5, out_type="gif", limit=None):
    lefts  = load_images_from_dir(left_dir)
    mids   = load_images_from_dir(mid_dir)
    rights = load_images_from_dir(right_dir)

    n = min(len(lefts), len(mids), len(rights))
    if limit:
        n = min(n, limit)

    Path(out_dir).mkdir(parents=True, exist_ok=True)
    out_file = os.path.join(out_dir, f"concat.{out_type}")

    frames = []
    for i in range(n):
        L = Image.open(lefts[i])
        M = Image.open(mids[i])
        R = Image.open(rights[i])
        frame = concat_h([L, M, R])
        frames.append(frame)

    if out_type == "gif":
        duration = 1.0 / max(fps, 0.0001)
        imageio.mimsave(out_file, [np.array(f) for f in frames], duration=duration, loop=0)
    elif out_type == "mp4":
        with imageio.get_writer(out_file, fps=fps, codec="libx264") as w:
            for f in frames:
                arr = np.array(f)
                arr = ensure_even_dims(arr)
                w.append_data(arr)
    else:
        raise ValueError("지원되는 출력 형식은 gif 또는 mp4")

    print(f"저장 완료: {out_file}  (fps={fps}, frames={len(frames)})")

# -------- 실행 예시 --------
if __name__ == "__main__":
    # ===== 사용자 입력 =====
    left_dir  = "C:/Users/ust21/shjo/MCC/cmems_fig/chl/2010-2013/"
    mid_dir   = "D:/shjo/nifs02_clm/fennel_2010/figs/chl/"
    right_dir = "D:/shjo/nifs02_clm/npzd_2010/figs/chl/"
    out_dir   = "C:/Users/ust21/shjo/MCC/compare_clm/"
    fps       = 2           # 프레임 속도
    out_type  = "mp4"       # "gif" 또는 "mp4"
    limit     = None        # 프레임 수 제한 (없으면 전체)

    make_animation(left_dir, mid_dir, right_dir, out_dir, fps=fps, out_type=out_type, limit=limit)
