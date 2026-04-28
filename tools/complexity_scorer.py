"""
合成难度评分器 - Complexity Scorer
综合评估序列的合成难度（0-100分）
"""

from tools.sequence_validator import calc_gc_content, basic_stats
from tools.structure_analyzer import find_g_quadruplex, find_z_dna, find_i_motif, find_h_dna


def score_complexity(seq: str) -> dict:
    """
    综合合成难度评分 (0-100)。
    考虑因子：GC极端区、重复序列、特殊结构密度、同源区段、序列长度。

    Args:
        seq: DNA 序列
    Returns:
        dict: {score, grade, factors, difficult_regions}
    """
    seq = seq.upper().strip()
    n = len(seq)
    if n == 0:
        return {"error": "Empty sequence"}

    factors = {}
    total_penalty = 0

    # 1. 长度因子 (0-15分)
    if n > 100000:
        length_penalty = 15
    elif n > 50000:
        length_penalty = 12
    elif n > 10000:
        length_penalty = 8
    elif n > 5000:
        length_penalty = 5
    else:
        length_penalty = 2
    factors["length"] = {"penalty": length_penalty, "value": n, "note": f"Sequence length: {n}bp"}
    total_penalty += length_penalty

    # 2. GC 含量因子 (0-20分)
    gc_data = calc_gc_content(seq, window=100)
    global_gc = gc_data["global_gc"]
    gc_penalty = 0
    if global_gc > 0.65 or global_gc < 0.30:
        gc_penalty = 15
    elif global_gc > 0.60 or global_gc < 0.35:
        gc_penalty = 8
    elif global_gc > 0.55 or global_gc < 0.40:
        gc_penalty = 3

    # 极端 GC 区段加分
    extreme_regions = len(gc_data["gc_regions"].get("high_gc", [])) + len(gc_data["gc_regions"].get("low_gc", []))
    gc_penalty += min(extreme_regions * 2, 5)
    gc_penalty = min(gc_penalty, 20)
    factors["gc_content"] = {"penalty": gc_penalty, "value": round(global_gc, 3),
                             "extreme_regions": extreme_regions}
    total_penalty += gc_penalty

    # 3. 重复序列因子 (0-20分)
    stats = basic_stats(seq)
    repeats = stats.get("simple_repeats", [])
    long_repeats = [r for r in repeats if r["total_length"] >= 20]
    repeat_penalty = min(len(long_repeats) * 3, 20)
    factors["repeats"] = {"penalty": repeat_penalty, "long_repeats": len(long_repeats),
                          "total_repeats": len(repeats)}
    total_penalty += repeat_penalty

    # 4. G-quadruplex 因子 (0-15分)
    g4 = find_g_quadruplex(seq)
    g4_high = sum(1 for s in g4.get("g4_sites", []) if s.get("stability") == "high")
    g4_penalty = min(g4_high * 5 + (g4["g4_count"] - g4_high) * 2, 15)
    factors["g_quadruplex"] = {"penalty": g4_penalty, "count": g4["g4_count"],
                               "high_stability": g4_high}
    total_penalty += g4_penalty

    # 5. Z-DNA 因子 (0-10分)
    zdna = find_z_dna(seq)
    zdna_high = sum(1 for r in zdna.get("z_dna_regions", []) if r.get("z_propensity") == "high")
    zdna_penalty = min(zdna_high * 3 + (zdna["z_dna_count"] - zdna_high), 10)
    factors["z_dna"] = {"penalty": zdna_penalty, "count": zdna["z_dna_count"]}
    total_penalty += zdna_penalty

    # 6. i-Motif / H-DNA 因子 (0-10分)
    imotif = find_i_motif(seq)
    hdna = find_h_dna(seq)
    other_penalty = min(imotif["i_motif_count"] * 2 + hdna["h_dna_count"], 10)
    factors["other_structures"] = {"penalty": other_penalty,
                                   "i_motif": imotif["i_motif_count"],
                                   "h_dna": hdna["h_dna_count"]}
    total_penalty += other_penalty

    # 7. 同源区段检测 (0-10分) - 检查序列内部自同源
    homo_penalty = 0
    check_len = min(n, 20000)  # 限制计算量
    window = 50
    for i in range(0, check_len - window, window):
        seg = seq[i:i + window]
        # 在后面的序列中搜索
        rest = seq[i + window + 100:]
        if seg in rest:
            homo_penalty += 2
    homo_penalty = min(homo_penalty, 10)
    factors["homology"] = {"penalty": homo_penalty, "note": "Internal sequence homology"}
    total_penalty += homo_penalty

    score = min(total_penalty, 100)

    if score >= 70:
        grade = "EXTREMELY_DIFFICULT"
    elif score >= 50:
        grade = "DIFFICULT"
    elif score >= 30:
        grade = "MODERATE"
    elif score >= 15:
        grade = "EASY"
    else:
        grade = "TRIVIAL"

    return {
        "complexity_score": score,
        "grade": grade,
        "factors": factors,
        "recommendation": _get_recommendation(grade, factors),
    }


def identify_difficult_regions(seq: str) -> dict:
    """标记具体困难区段及原因。"""
    seq = seq.upper().strip()
    regions = []

    # GC 极端区
    gc_data = calc_gc_content(seq, window=100)
    for r in gc_data["gc_regions"].get("high_gc", []):
        regions.append({"start": r["start"], "end": r["end"], "type": "high_gc",
                       "detail": f"GC={r['gc']}", "severity": "HIGH"})
    for r in gc_data["gc_regions"].get("low_gc", []):
        regions.append({"start": r["start"], "end": r["end"], "type": "low_gc",
                       "detail": f"GC={r['gc']}", "severity": "MODERATE"})

    # G4
    g4 = find_g_quadruplex(seq)
    for s in g4.get("g4_sites", []):
        regions.append({"start": s["position"], "end": s["end"], "type": "G-quadruplex",
                       "detail": s.get("stability", ""), "severity": s.get("synthesis_risk", "MODERATE")})

    # Z-DNA
    zdna = find_z_dna(seq)
    for r in zdna.get("z_dna_regions", []):
        regions.append({"start": r["position"], "end": r["end"], "type": "Z-DNA",
                       "detail": f"propensity={r['z_propensity']}", "severity": "MODERATE"})

    # 排序
    regions.sort(key=lambda r: r["start"])

    return {
        "total_difficult_regions": len(regions),
        "regions": regions,
        "coverage": _calc_coverage(regions, len(seq)),
    }


def _calc_coverage(regions, seq_len):
    if seq_len == 0:
        return 0
    covered = set()
    for r in regions:
        for i in range(r["start"], min(r["end"], seq_len)):
            covered.add(i)
    return round(len(covered) / seq_len, 4)


def _get_recommendation(grade, factors):
    recs = []
    if grade in ("EXTREMELY_DIFFICULT", "DIFFICULT"):
        recs.append("Consider splitting into smaller fragments with optimized overlaps.")
    if factors.get("gc_content", {}).get("penalty", 0) > 10:
        recs.append("GC content is extreme. Use codon optimization or non-coding region modification.")
    if factors.get("g_quadruplex", {}).get("high_stability", 0) > 0:
        recs.append("G-quadruplex detected. Add 5% DMSO or 1M Betaine to PCR. Use high-fidelity polymerase.")
    if factors.get("repeats", {}).get("long_repeats", 0) > 3:
        recs.append("Long repeats detected. Consider breaking repeats with synonymous substitutions.")
    if not recs:
        recs.append("Sequence is within normal synthesis parameters.")
    return recs
