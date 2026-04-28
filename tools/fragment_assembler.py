"""
片段组装器 - Fragment Assembler
Gibson Assembly / Golden Gate 片段拆分、Overlap 优化、组装报告
"""

import json
import os
from typing import Optional
from tools.thermodynamics import calc_tm, calc_nn_thermodynamics

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")


def split_gibson(
    seq: str,
    target_size: int = 5000,
    overlap: int = 40,
    na: float = 50.0,
    mg: float = 1.5,
) -> dict:
    """
    Gibson Assembly 策略拆分序列。

    Args:
        seq: DNA 序列
        target_size: 目标片段大小 (bp)
        overlap: 重叠区长度 (bp)
        na: [Na+] mM
        mg: [Mg2+] mM
    """
    seq = seq.upper().strip()
    n = len(seq)

    if n < target_size:
        return {
            "fragments": [{"name": "F1", "start": 0, "end": n, "length": n, "sequence": seq}],
            "total_fragments": 1,
            "note": "Sequence smaller than target size, no splitting needed.",
        }

    fragments = []
    pos = 0
    frag_id = 1

    while pos < n:
        end = min(pos + target_size, n)

        # 尝试在 ±200bp 范围内找最优切割点
        best_cut = end
        best_score = -1

        if end < n:
            for offset in range(-200, 201, 5):
                cut = end + offset
                if cut <= pos + 1000 or cut >= n:
                    continue

                # 评估切割点处的 overlap 质量
                ol_start = max(0, cut - overlap)
                ol_seq = seq[ol_start:cut]
                if len(ol_seq) < 20:
                    continue

                gc = sum(1 for b in ol_seq if b in "GC") / len(ol_seq)
                # 理想 GC: 40-60%
                gc_score = 1.0 - abs(gc - 0.50) * 4

                # 避免同碱基 run
                max_run = 1
                run = 1
                for i in range(1, len(ol_seq)):
                    if ol_seq[i] == ol_seq[i-1]:
                        run += 1
                        max_run = max(max_run, run)
                    else:
                        run = 1
                run_score = 1.0 if max_run <= 3 else 0.5 if max_run <= 5 else 0.0

                score = gc_score + run_score
                if score > best_score:
                    best_score = score
                    best_cut = cut

        frag_seq = seq[pos:best_cut]

        # 计算 overlap
        overlap_5 = ""
        overlap_3 = ""
        if frag_id > 1:
            overlap_5 = seq[max(0, pos - overlap):pos]
        if best_cut < n:
            overlap_3 = seq[best_cut - overlap:best_cut]

        # Tm 计算
        tm_5 = calc_tm(overlap_5, na, mg).get("Tm_C", 0) if overlap_5 and len(overlap_5) >= 10 else None
        tm_3 = calc_tm(overlap_3, na, mg).get("Tm_C", 0) if overlap_3 and len(overlap_3) >= 10 else None

        fragments.append({
            "name": f"F{frag_id}",
            "start": pos,
            "end": best_cut,
            "length": best_cut - pos,
            "gc_content": round(sum(1 for b in frag_seq if b in "GC") / max(len(frag_seq), 1), 3),
            "overlap_5prime": overlap_5,
            "overlap_3prime": overlap_3,
            "overlap_5_Tm_C": tm_5,
            "overlap_3_Tm_C": tm_3,
        })

        pos = best_cut
        frag_id += 1

    return {
        "method": "Gibson Assembly",
        "total_fragments": len(fragments),
        "target_fragment_size": target_size,
        "overlap_length": overlap,
        "fragments": fragments,
        "conditions": {"Na_mM": na, "Mg_mM": mg},
    }


def split_golden_gate(
    seq: str,
    target_size: int = 5000,
    enzyme: str = "BsaI",
) -> dict:
    """
    Golden Gate Assembly 策略拆分序列。

    Args:
        seq: DNA 序列
        target_size: 目标片段大小 (bp)
        enzyme: Type IIS 限制性酶
    """
    seq = seq.upper().strip()
    n = len(seq)

    with open(os.path.join(DATA_DIR, "restriction_enzymes.json"), "r") as f:
        enzyme_db = json.load(f)["enzymes"]

    if enzyme not in enzyme_db:
        return {"error": f"Unknown enzyme: {enzyme}"}

    enz_info = enzyme_db[enzyme]

    fragments = []
    pos = 0
    frag_id = 1
    used_overhangs = set()

    while pos < n:
        end = min(pos + target_size, n)
        if end >= n:
            end = n

        # 生成 4bp overhang
        if end < n:
            overhang = seq[end:end + 4] if end + 4 <= n else "NNNN"
            # 确保 overhang 唯一
            attempts = 0
            while overhang in used_overhangs and attempts < 20:
                end += 1
                overhang = seq[end:end + 4] if end + 4 <= n else "NNNN"
                attempts += 1
            used_overhangs.add(overhang)
        else:
            overhang = None

        fragments.append({
            "name": f"F{frag_id}",
            "start": pos,
            "end": end,
            "length": end - pos,
            "overhang_3prime": overhang,
            "gc_content": round(
                sum(1 for b in seq[pos:end] if b in "GC") / max(end - pos, 1), 3
            ),
        })

        pos = end
        frag_id += 1

    return {
        "method": "Golden Gate Assembly",
        "enzyme": enzyme,
        "enzyme_info": enz_info,
        "total_fragments": len(fragments),
        "fragments": fragments,
        "unique_overhangs": len(used_overhangs),
        "note": f"Ensure no internal {enzyme} sites exist in fragments.",
    }


def optimize_overlaps(fragments: list, target_tm: float = 62.0, na: float = 50.0, mg: float = 1.5) -> dict:
    """
    优化 Overlap Tm 均衡性。调整 overlap 长度使所有 Tm 接近目标值。

    Args:
        fragments: 片段列表 [{name, overlap_3prime, ...}]
        target_tm: 目标 Tm (°C)
    """
    optimized = []
    for frag in fragments:
        ol = frag.get("overlap_3prime", "")
        if not ol or len(ol) < 10:
            optimized.append({**frag, "optimized": False})
            continue

        current_tm = calc_tm(ol, na, mg).get("Tm_C", 0)
        best_ol = ol
        best_diff = abs(current_tm - target_tm)

        # 尝试增减 overlap 长度
        for delta in range(-10, 11, 2):
            new_len = len(ol) + delta
            if new_len < 15 or new_len > 60:
                continue
            # 这里需要原始序列，简化处理
            test_ol = ol[:new_len] if delta < 0 else ol
            tm = calc_tm(test_ol, na, mg).get("Tm_C", 0)
            diff = abs(tm - target_tm)
            if diff < best_diff:
                best_diff = diff
                best_ol = test_ol

        optimized.append({
            **frag,
            "overlap_3prime": best_ol,
            "overlap_length_optimized": len(best_ol),
            "Tm_optimized_C": calc_tm(best_ol, na, mg).get("Tm_C", 0),
            "optimized": True,
        })

    return {"fragments": optimized, "target_tm_C": target_tm}


def generate_assembly_report(fragments: list, method: str = "Gibson") -> dict:
    """
    生成完整的组装方案报告。

    Args:
        fragments: 拆分后的片段列表
        method: 组装方法
    """
    report = {
        "method": method,
        "total_fragments": len(fragments),
        "fragment_summary": [],
        "protocol": {},
    }

    total_length = 0
    tms = []

    for frag in fragments:
        total_length += frag.get("length", 0)
        summary = {
            "name": frag.get("name", ""),
            "length": frag.get("length", 0),
            "gc_content": frag.get("gc_content", 0),
        }
        if method == "Gibson":
            ol_tm = frag.get("overlap_3_Tm_C") or frag.get("Tm_optimized_C")
            if ol_tm:
                tms.append(ol_tm)
                summary["overlap_Tm"] = ol_tm
        report["fragment_summary"].append(summary)

    if method == "Gibson":
        report["protocol"] = {
            "master_mix": "Gibson Assembly Master Mix (NEB E2611)",
            "incubation": "50°C for 60 min",
            "fragment_ratio": "equimolar, 50-100 ng each for 2-3 fragments",
            "total_volume_uL": 20,
            "competent_cells": "NEB 10-beta or DH5α",
        }
        if tms:
            import statistics
            report["overlap_tm_stats"] = {
                "mean_C": round(statistics.mean(tms), 1),
                "std_C": round(statistics.stdev(tms) if len(tms) > 1 else 0, 1),
                "range_C": round(max(tms) - min(tms), 1),
            }
    elif method == "Golden Gate":
        report["protocol"] = {
            "reaction": "30 cycles of (37°C 5min, 16°C 5min), then 55°C 5min, 80°C 5min",
            "enzyme_amount": "1 µL Type IIS RE + 1 µL T4 DNA Ligase per 20 µL reaction",
            "fragment_amount": "75 ng each fragment",
        }

    report["total_assembled_length"] = total_length
    return report


def suggest_pcr_protocol(fragment: dict) -> dict:
    """
    为特定片段建议 PCR 条件。

    Args:
        fragment: 片段信息 {length, gc_content, ...}
    """
    length = fragment.get("length", 0)
    gc = fragment.get("gc_content", 0.50)

    protocol = {
        "polymerase": "Q5 High-Fidelity (NEB M0491)",
        "initial_denaturation": {"temp_C": 98, "time_s": 30},
        "cycles": 30,
        "denaturation": {"temp_C": 98, "time_s": 10},
        "annealing": {"temp_C": 65, "time_s": 20},
        "extension": {"temp_C": 72, "time_s": max(30, length // 1000 * 30)},
        "final_extension": {"temp_C": 72, "time_s": 120},
        "additives": [],
    }

    # GC 调整
    if gc > 0.65:
        protocol["additives"].append({"name": "DMSO", "concentration_pct": 5})
        protocol["initial_denaturation"]["time_s"] = 60
        protocol["annealing"]["temp_C"] = 68
    elif gc > 0.60:
        protocol["additives"].append({"name": "Betaine", "concentration_M": 1.0})

    # 低 GC 调整
    if gc < 0.30:
        protocol["annealing"]["temp_C"] = 58
        protocol["extension"]["time_s"] = int(protocol["extension"]["time_s"] * 1.5)

    # 长片段调整
    if length > 5000:
        protocol["polymerase"] = "Q5 High-Fidelity or PrimeSTAR GXL (for >5kb)"
        protocol["extension"]["time_s"] = max(60, length // 1000 * 60)
        protocol["cycles"] = 25

    if length > 10000:
        protocol["note"] = "Consider using long-range PCR kit. Fragment >10kb may require optimization."

    return {"fragment": fragment.get("name", ""), "protocol": protocol}
