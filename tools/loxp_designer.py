"""
loxP 位点设计器 - LoxP Designer
loxP/lox511/lox2272 位点布局、正交性验证、SCRaMbLE 模拟
"""

import random
from typing import Optional

# loxP 及其变体序列
LOXP_VARIANTS = {
    "loxP": {
        "sequence": "ATAACTTCGTATAATGTATGCTATACGAAGTTAT",
        "left_arm": "ATAACTTCGTATA",
        "spacer": "ATGTATGC",
        "right_arm": "TATACGAAGTTAT",
        "description": "Wild-type loxP site (34bp)",
    },
    "lox511": {
        "sequence": "ATAACTTCGTATAAAGTATGCTATACGAAGTTAT",
        "left_arm": "ATAACTTCGTATA",
        "spacer": "AAGTATGC",
        "right_arm": "TATACGAAGTTAT",
        "description": "lox511 - orthogonal to loxP, single bp change in spacer",
    },
    "lox2272": {
        "sequence": "ATAACTTCGTATAAAGTATACTATACGAAGTTAT",
        "left_arm": "ATAACTTCGTATA",
        "spacer": "AAGTATAC",
        "right_arm": "TATACGAAGTTAT",
        "description": "lox2272 - orthogonal to both loxP and lox511",
    },
    "loxFAS": {
        "sequence": "ATAACTTCGTATATTCTATACTATACGAAGTTAT",
        "spacer": "TTCTATAC",
        "description": "loxFAS - asymmetric spacer, irreversible excision",
    },
    "lox5171": {
        "sequence": "ATAACTTCGTATAATGTGTACTATACGAAGTTAT",
        "spacer": "ATGTGTAC",
        "description": "lox5171 - orthogonal variant",
    },
}


def design_loxp_sites(
    seq: str,
    interval: int = 5000,
    variant: str = "loxP",
    avoid_regions: Optional[list] = None,
    min_distance_from_feature: int = 100,
) -> dict:
    """
    在序列中自动布置 loxP 位点。

    Args:
        seq: DNA 序列
        interval: 间隔 (bp)，默认每5kb一个
        variant: loxP 变体类型
        avoid_regions: 需要避开的区域 [{start, end, name}]（如必需基因）
        min_distance_from_feature: 与 feature 的最小距离
    """
    seq = seq.upper().strip()
    n = len(seq)
    avoid = avoid_regions or []

    if variant not in LOXP_VARIANTS:
        return {"error": f"Unknown variant: {variant}. Available: {list(LOXP_VARIANTS.keys())}"}

    loxp_seq = LOXP_VARIANTS[variant]["sequence"]
    sites = []
    target_positions = list(range(interval, n - 100, interval))

    for target in target_positions:
        # 在 target ± 500bp 范围内寻找最优位点
        best_pos = target
        best_score = -1

        for offset in range(-500, 501, 10):
            pos = target + offset
            if pos < 50 or pos + len(loxp_seq) > n - 50:
                continue

            # 检查是否在 avoid 区域内
            in_avoid = False
            for ar in avoid:
                if ar["start"] - min_distance_from_feature <= pos <= ar["end"] + min_distance_from_feature:
                    in_avoid = True
                    break
            if in_avoid:
                continue

            # 评分：离 avoid 区域越远越好，离目标越近越好
            dist_to_avoid = min(
                (abs(pos - ar["start"]) + abs(pos - ar["end"])) for ar in avoid
            ) if avoid else 10000

            dist_to_target = abs(pos - target)
            score = dist_to_avoid - dist_to_target * 0.5

            # 检查局部 GC 含量（避免极端区域）
            local_seq = seq[max(0, pos - 50):pos + len(loxp_seq) + 50]
            local_gc = sum(1 for b in local_seq if b in "GC") / len(local_seq)
            if 0.35 < local_gc < 0.65:
                score += 100

            if score > best_score:
                best_score = score
                best_pos = pos

        sites.append({
            "position": best_pos,
            "variant": variant,
            "sequence": loxp_seq,
            "local_gc": round(
                sum(1 for b in seq[max(0, best_pos-50):best_pos+len(loxp_seq)+50] if b in "GC")
                / min(len(loxp_seq) + 100, n), 3
            ),
            "distance_to_nearest_avoid": min(
                min(abs(best_pos - ar["start"]), abs(best_pos - ar["end"])) for ar in avoid
            ) if avoid else None,
        })

    return {
        "variant": variant,
        "variant_info": LOXP_VARIANTS[variant],
        "total_sites": len(sites),
        "interval_bp": interval,
        "sites": sites,
        "modified_sequence_length": n + len(sites) * len(loxp_seq),
    }


def check_orthogonality(variants: list = None) -> dict:
    """
    验证 loxP 变体间的正交性。
    正交性由 spacer 序列决定 - 不同 spacer 的 lox 位点不会交叉重组。

    Args:
        variants: 要检查的变体列表，默认全部
    """
    if variants is None:
        variants = list(LOXP_VARIANTS.keys())

    pairs = []
    for i, v1 in enumerate(variants):
        if v1 not in LOXP_VARIANTS:
            continue
        for v2 in variants[i + 1:]:
            if v2 not in LOXP_VARIANTS:
                continue
            s1 = LOXP_VARIANTS[v1].get("spacer", "")
            s2 = LOXP_VARIANTS[v2].get("spacer", "")

            # 比较 spacer 序列
            mismatches = sum(1 for a, b in zip(s1, s2) if a != b) if len(s1) == len(s2) else len(s1)
            orthogonal = mismatches >= 1  # 至少1个碱基不同
            recombination_risk = "NONE" if mismatches >= 3 else "LOW" if mismatches >= 2 else "MODERATE" if mismatches >= 1 else "HIGH"

            pairs.append({
                "variant_1": v1, "variant_2": v2,
                "spacer_1": s1, "spacer_2": s2,
                "mismatches": mismatches,
                "orthogonal": orthogonal,
                "recombination_risk": recombination_risk,
            })

    all_orthogonal = all(p["orthogonal"] for p in pairs)
    return {
        "variants_checked": variants,
        "pairs": pairs,
        "all_orthogonal": all_orthogonal,
        "recommendation": "All variants are orthogonal and safe to use together."
            if all_orthogonal else "WARNING: Some pairs may cross-recombine.",
    }


def simulate_scramble(
    seq_length: int,
    loxp_positions: list,
    n_events: int = 100,
    seed: int = 42,
) -> dict:
    """
    模拟 Cre 介导的 SCRaMbLE 拓扑演变。

    Args:
        seq_length: 序列总长度
        loxp_positions: loxP 位点位置列表
        n_events: 模拟的重组事件数
        seed: 随机种子
    """
    random.seed(seed)
    positions = sorted(loxp_positions)
    n_sites = len(positions)

    if n_sites < 2:
        return {"error": "Need at least 2 loxP sites for SCRaMbLE simulation"}

    # 定义 segments（loxP 之间的区段）
    segments = []
    for i in range(n_sites - 1):
        segments.append({
            "id": i,
            "start": positions[i],
            "end": positions[i + 1],
            "length": positions[i + 1] - positions[i],
            "orientation": "+",
        })

    events = []
    current_segments = list(range(len(segments)))
    deletions = 0
    inversions = 0

    for event_id in range(n_events):
        if len(current_segments) < 2:
            break

        # 随机选两个 loxP 位点
        idx1, idx2 = sorted(random.sample(range(len(current_segments)), 2))

        # 相同方向 → 删除；反向 → 反转
        event_type = random.choice(["deletion", "inversion"])

        if event_type == "deletion":
            deleted = current_segments[idx1:idx2]
            current_segments = current_segments[:idx1] + current_segments[idx2:]
            deletions += 1
            events.append({
                "event_id": event_id,
                "type": "deletion",
                "segments_affected": deleted,
                "remaining_segments": len(current_segments),
            })
        else:
            inverted = current_segments[idx1:idx2][::-1]
            current_segments = current_segments[:idx1] + inverted + current_segments[idx2:]
            inversions += 1
            events.append({
                "event_id": event_id,
                "type": "inversion",
                "segments_affected": list(range(idx1, idx2)),
                "remaining_segments": len(current_segments),
            })

    # 计算最终序列长度
    remaining_length = sum(segments[i]["length"] for i in current_segments if i < len(segments))

    return {
        "initial_segments": len(segments),
        "initial_loxp_sites": n_sites,
        "total_events": len(events),
        "deletions": deletions,
        "inversions": inversions,
        "final_segments": len(current_segments),
        "final_length_estimate": remaining_length,
        "length_reduction_pct": round((1 - remaining_length / max(seq_length, 1)) * 100, 1),
        "events_log": events[:20],  # 只返回前20个事件
        "topology_summary": f"{deletions} deletions, {inversions} inversions over {len(events)} events",
    }
