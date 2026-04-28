"""
序列验证工具 - Sequence Validator
GC% 计算、长度检查、限制性酶切位点扫描
"""

import json
import re
import os
from typing import Optional

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

IUPAC_DNA = set("ATCGRYSWKMBDHVNatcgryswkmdbhvn")


def validate_sequence(seq: str) -> dict:
    """
    验证 DNA 序列的合法性。
    
    Args:
        seq: DNA 序列字符串
    
    Returns:
        dict: {valid: bool, length: int, invalid_chars: list, invalid_positions: list}
    """
    seq = seq.strip().replace("\n", "").replace("\r", "").replace(" ", "")
    invalid_chars = []
    invalid_positions = []
    for i, c in enumerate(seq):
        if c not in IUPAC_DNA:
            invalid_chars.append(c)
            invalid_positions.append(i)
    return {
        "valid": len(invalid_chars) == 0,
        "length": len(seq),
        "invalid_chars": list(set(invalid_chars)),
        "invalid_positions": invalid_positions[:20],  # 最多报20个
        "cleaned_sequence_length": len(seq),
    }


def calc_gc_content(seq: str, window: int = 0) -> dict:
    """
    计算 GC 含量（全局 & 滑窗）。
    
    Args:
        seq: DNA 序列
        window: 滑窗大小，0 表示仅计算全局
    
    Returns:
        dict: {global_gc: float, window_gc: list[{pos, gc}], gc_regions: dict}
    """
    seq = seq.upper().replace("\n", "").replace(" ", "")
    total = len(seq)
    if total == 0:
        return {"global_gc": 0, "window_gc": [], "gc_regions": {}}

    gc_count = seq.count("G") + seq.count("C")
    global_gc = gc_count / total

    window_gc = []
    if window > 0 and total >= window:
        for i in range(0, total - window + 1, max(1, window // 10)):
            w = seq[i : i + window]
            wgc = (w.count("G") + w.count("C")) / window
            window_gc.append({"position": i, "gc": round(wgc, 4)})

    # 识别极端 GC 区段
    high_gc_regions = []
    low_gc_regions = []
    scan_window = min(100, total)
    if total >= scan_window:
        for i in range(0, total - scan_window + 1, scan_window // 2):
            w = seq[i : i + scan_window]
            wgc = (w.count("G") + w.count("C")) / scan_window
            if wgc > 0.70:
                high_gc_regions.append({"start": i, "end": i + scan_window, "gc": round(wgc, 4)})
            elif wgc < 0.25:
                low_gc_regions.append({"start": i, "end": i + scan_window, "gc": round(wgc, 4)})

    return {
        "global_gc": round(global_gc, 4),
        "window_gc": window_gc,
        "gc_regions": {
            "high_gc": high_gc_regions,
            "low_gc": low_gc_regions,
        },
    }


def find_restriction_sites(seq: str, enzymes: Optional[list] = None) -> dict:
    """
    扫描限制性酶切位点。
    
    Args:
        seq: DNA 序列
        enzymes: 指定酶列表，None 表示扫描全部
    
    Returns:
        dict: {enzyme_name: [{position, sequence}]}
    """
    seq = seq.upper().replace("\n", "").replace(" ", "")
    
    with open(os.path.join(DATA_DIR, "restriction_enzymes.json"), "r") as f:
        enzyme_data = json.load(f)

    results = {}
    enzyme_db = enzyme_data["enzymes"]
    target_enzymes = enzymes if enzymes else list(enzyme_db.keys())

    iupac_map = {
        "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
        "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
        "H": "[ACT]", "V": "[ACG]", "N": "[ATCG]",
    }

    for enz_name in target_enzymes:
        if enz_name not in enzyme_db:
            continue
        rec_seq = enzyme_db[enz_name]["sequence"]
        # 将 IUPAC 简并碱基转为正则
        pattern = ""
        for c in rec_seq:
            pattern += iupac_map.get(c, c)

        sites = []
        for m in re.finditer(pattern, seq):
            sites.append({
                "position": m.start(),
                "sequence": m.group(),
            })
        if sites:
            results[enz_name] = {
                "recognition_sequence": rec_seq,
                "sites": sites,
                "count": len(sites),
            }

    return {
        "total_enzymes_scanned": len(target_enzymes),
        "enzymes_with_sites": len(results),
        "results": results,
    }


def basic_stats(seq: str) -> dict:
    """
    序列基本统计信息。
    
    Args:
        seq: DNA 序列
    
    Returns:
        dict: 综合统计
    """
    seq = seq.upper().replace("\n", "").replace(" ", "")
    total = len(seq)
    if total == 0:
        return {"error": "Empty sequence"}

    counts = {base: seq.count(base) for base in "ATCG"}
    gc = counts["G"] + counts["C"]
    at = counts["A"] + counts["T"]

    # AT/GC skew
    at_skew = (counts["A"] - counts["T"]) / max(at, 1)
    gc_skew = (counts["G"] - counts["C"]) / max(gc, 1)

    # 简单重复序列扫描
    repeats = []
    for unit_len in range(1, 7):
        for i in range(total - unit_len * 3):
            unit = seq[i : i + unit_len]
            count = 1
            j = i + unit_len
            while j + unit_len <= total and seq[j : j + unit_len] == unit:
                count += 1
                j += unit_len
            if count >= 4 and unit_len * count >= 12:
                repeats.append({
                    "position": i,
                    "unit": unit,
                    "count": count,
                    "total_length": unit_len * count,
                })

    # 去重
    seen = set()
    unique_repeats = []
    for r in repeats:
        key = (r["position"], r["unit"])
        if key not in seen:
            seen.add(key)
            unique_repeats.append(r)

    return {
        "length": total,
        "base_composition": counts,
        "gc_content": round(gc / total, 4),
        "at_skew": round(at_skew, 4),
        "gc_skew": round(gc_skew, 4),
        "simple_repeats": unique_repeats[:50],  # 最多报50个
        "size_category": (
            "small" if total < 1000
            else "medium" if total < 10000
            else "large" if total < 50000
            else "mega"
        ),
    }
