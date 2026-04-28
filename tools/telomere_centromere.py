"""
端粒与着丝粒设计器 - Telomere & Centromere Designer
端粒重复序列生成、着丝粒设计、稳定性评估
"""

from typing import Optional


# 各物种端粒重复序列
TELOMERE_REPEATS = {
    "yeast": {"unit": "TGTGGGTGTGGTG", "alt_unit": "TG1-3", "description": "S. cerevisiae telomere (irregular TG1-3 repeats)"},
    "human": {"unit": "TTAGGG", "description": "H. sapiens telomere"},
    "mouse": {"unit": "TTAGGG", "description": "M. musculus telomere (same as human)"},
    "arabidopsis": {"unit": "TTTAGGG", "description": "A. thaliana telomere"},
    "tetrahymena": {"unit": "TTGGGG", "description": "T. thermophila telomere"},
    "drosophila": {"unit": "HeT-A/TART", "description": "D. melanogaster uses retrotransposons, not simple repeats"},
}

# 酵母着丝粒元件 (S. cerevisiae point centromere)
YEAST_CENTROMERE = {
    "CDE_I": {
        "consensus": "RTCACRTG",  # R = purine
        "example": "ATCACATG",
        "length": 8,
        "description": "Centromere DNA Element I (8bp conserved)",
    },
    "CDE_II": {
        "description": "Centromere DNA Element II (78-86bp, >90% AT)",
        "at_content_min": 0.90,
        "length_range": [78, 86],
    },
    "CDE_III": {
        "consensus": "TGTTTXTGNTTTCCGAAANNNNNAAAA",
        "length": 26,
        "description": "Centromere DNA Element III (26bp, most critical for function)",
    },
}


def design_telomere(organism: str = "yeast", length: int = 2000, add_subtelomeric: bool = True) -> dict:
    """
    生成端粒重复序列。

    Args:
        organism: 目标物种
        length: 目标长度 (bp)
        add_subtelomeric: 是否添加亚端粒缓冲区
    """
    if organism not in TELOMERE_REPEATS:
        return {"error": f"Unknown organism. Available: {list(TELOMERE_REPEATS.keys())}"}

    info = TELOMERE_REPEATS[organism]

    if organism == "drosophila":
        return {"error": "Drosophila uses retrotransposons for telomeres, not simple repeats. Manual design required."}

    if organism == "yeast":
        # S. cerevisiae has irregular TG1-3 repeats
        import random
        random.seed(42)
        telomere_seq = ""
        while len(telomere_seq) < length:
            g_count = random.choice([1, 2, 3])
            telomere_seq += "T" + "G" * g_count
        telomere_seq = telomere_seq[:length]
    else:
        unit = info["unit"]
        repeats_needed = (length // len(unit)) + 1
        telomere_seq = (unit * repeats_needed)[:length]

    # 亚端粒缓冲区 (Y' element analog, ~300bp)
    subtelomeric = ""
    if add_subtelomeric and organism == "yeast":
        # 简化版 Y' element 核心区
        subtelomeric = "ATCGATCGATCG" * 25  # 300bp placeholder
        subtelomeric = subtelomeric[:300]

    # G4 风险评估
    g4_risk = "HIGH" if organism in ["human", "mouse", "tetrahymena"] else "MODERATE"

    result = {
        "organism": organism,
        "telomere_info": info,
        "telomere_sequence": telomere_seq,
        "telomere_length": len(telomere_seq),
        "gc_content": round(sum(1 for b in telomere_seq if b in "GC") / len(telomere_seq), 3),
        "g4_risk": g4_risk,
        "synthesis_warnings": [],
    }

    if g4_risk == "HIGH":
        result["synthesis_warnings"].append(
            "G-rich telomere repeats are prone to G-quadruplex formation. "
            "Recommend: 5% DMSO, reduce extension temperature, use specialized polymerase."
        )

    if add_subtelomeric and subtelomeric:
        result["subtelomeric_buffer"] = subtelomeric
        result["subtelomeric_length"] = len(subtelomeric)
        result["total_length"] = len(telomere_seq) + len(subtelomeric)

    return result


def design_centromere(organism: str = "yeast", chromosome_id: str = "synIII") -> dict:
    """
    设计着丝粒序列。

    Args:
        organism: 目标物种
        chromosome_id: 染色体编号
    """
    if organism != "yeast":
        return {
            "error": f"Currently only S. cerevisiae point centromeres are supported. "
                     f"Regional centromeres (human, etc.) require specialized design.",
            "note": "For regional centromeres, consider using alpha-satellite DNA arrays.",
        }

    import random
    random.seed(hash(chromosome_id) % 2**32)

    # CDE I
    cde_i = YEAST_CENTROMERE["CDE_I"]["example"]

    # CDE II: 78-86bp, >90% AT
    cde_ii_len = random.randint(78, 86)
    cde_ii = ""
    for _ in range(cde_ii_len):
        if random.random() < 0.92:  # ~92% AT
            cde_ii += random.choice("AT")
        else:
            cde_ii += random.choice("GC")

    # CDE III (26bp conserved)
    cde_iii_template = "TGTTTTTGCTTTCCGAAA"
    cde_iii_variable = "".join(random.choice("ATCG") for _ in range(5))
    cde_iii = cde_iii_template + cde_iii_variable + "AAAA"

    # 组装
    centromere_seq = cde_i + cde_ii + cde_iii
    at_content = sum(1 for b in centromere_seq if b in "AT") / len(centromere_seq)

    return {
        "organism": organism,
        "chromosome_id": chromosome_id,
        "centromere_sequence": centromere_seq,
        "centromere_length": len(centromere_seq),
        "elements": {
            "CDE_I": {"sequence": cde_i, "length": len(cde_i)},
            "CDE_II": {"sequence": cde_ii, "length": len(cde_ii), "at_content": round(
                sum(1 for b in cde_ii if b in "AT") / len(cde_ii), 3)},
            "CDE_III": {"sequence": cde_iii, "length": len(cde_iii)},
        },
        "overall_at_content": round(at_content, 3),
        "synthesis_warnings": [
            "High AT content (>90%) in CDE II may cause polymerase slippage.",
            "Consider flanking with GC-clamp sequences for PCR amplification.",
        ] if at_content > 0.85 else [],
    }


def assess_stability(element_seq: str, host: str = "yeast") -> dict:
    """
    评估功能元件在宿主中的稳定性。

    Args:
        element_seq: 元件 DNA 序列
        host: 宿主物种
    """
    seq = element_seq.upper().strip()
    n = len(seq)

    # 基本统计
    gc = sum(1 for b in seq if b in "GC") / max(n, 1)
    at = 1 - gc

    # 重复序列密度
    repeat_score = 0
    for unit_len in range(1, 7):
        for i in range(n - unit_len * 3):
            unit = seq[i:i + unit_len]
            count = 1
            j = i + unit_len
            while j + unit_len <= n and seq[j:j + unit_len] == unit:
                count += 1
                j += unit_len
            if count >= 4:
                repeat_score += count * unit_len
            if repeat_score > 1000:
                break
        if repeat_score > 1000:
            break

    repeat_density = min(repeat_score / max(n, 1), 1.0)

    # 稳定性评分
    stability_score = 100

    # GC 偏差惩罚
    if gc < 0.25 or gc > 0.75:
        stability_score -= 30
    elif gc < 0.35 or gc > 0.65:
        stability_score -= 15

    # 重复密度惩罚
    stability_score -= int(repeat_density * 40)

    # 长度惩罚（太长的重复元件不稳定）
    if n > 5000:
        stability_score -= 10
    if n > 10000:
        stability_score -= 10

    stability_score = max(0, stability_score)

    grade = ("STABLE" if stability_score >= 70
             else "MODERATE" if stability_score >= 40
             else "UNSTABLE")

    return {
        "element_length": n,
        "host": host,
        "gc_content": round(gc, 3),
        "repeat_density": round(repeat_density, 3),
        "stability_score": stability_score,
        "stability_grade": grade,
        "recommendations": _stability_recommendations(gc, repeat_density, grade),
    }


def optimize_for_synthesis(seq: str) -> dict:
    """
    针对合成困难区的优化建议。

    Args:
        seq: DNA 序列
    """
    seq = seq.upper().strip()
    suggestions = []

    gc = sum(1 for b in seq if b in "GC") / max(len(seq), 1)

    if gc > 0.70:
        suggestions.append({
            "type": "high_gc",
            "action": "Add 5-10% DMSO or 1M Betaine to PCR reactions",
            "protocol": {"DMSO_pct": 7, "denaturation_temp_C": 98, "denaturation_time_s": 30},
        })
    if gc < 0.25:
        suggestions.append({
            "type": "low_gc",
            "action": "Reduce annealing temperature by 3-5°C, extend extension time",
            "protocol": {"annealing_offset_C": -4, "extension_time_multiplier": 1.5},
        })

    # 检查连续同碱基
    for base in "ATCG":
        for i in range(len(seq) - 5):
            if seq[i:i+6] == base * 6:
                suggestions.append({
                    "type": "homopolymer",
                    "position": i,
                    "base": base,
                    "action": f"Homopolymer run of {base}×6+ at position {i}. Consider synonymous substitution if in CDS.",
                })
                break

    if not suggestions:
        suggestions.append({"type": "none", "action": "Sequence is within normal synthesis parameters."})

    return {"sequence_length": len(seq), "gc_content": round(gc, 3), "optimization_suggestions": suggestions}


def _stability_recommendations(gc, repeat_density, grade):
    recs = []
    if grade == "UNSTABLE":
        recs.append("Consider codon-optimizing repetitive regions to reduce internal homology.")
        recs.append("Add unique flanking sequences to anchor the element during replication.")
    if gc < 0.30:
        recs.append("Low GC may cause replication slippage. Intersperse GC-clamp sequences.")
    if repeat_density > 0.3:
        recs.append("High repeat density may cause recombination-mediated deletion in vivo.")
    if not recs:
        recs.append("Element appears stable for the target host.")
    return recs
