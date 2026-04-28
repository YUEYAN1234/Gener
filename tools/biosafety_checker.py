"""
生物安全检查器 - Biosafety Checker
管制病原体序列筛查、毒力因子扫描
"""

import json
import os
from typing import Optional

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")


class BiosafetyViolation(Exception):
    """生物安全红线触发异常"""
    def __init__(self, pathogen, risk_group, gene, identity):
        self.pathogen = pathogen
        self.risk_group = risk_group
        self.gene = gene
        self.identity = identity
        super().__init__(
            f"🚨 BIOSAFETY VIOLATION: Sequence matches {pathogen} ({gene}), "
            f"Risk Group {risk_group}, Identity {identity:.1%}"
        )


def _hamming_similarity(s1: str, s2: str) -> float:
    """计算两个等长序列的相似度"""
    if len(s1) != len(s2):
        return 0.0
    matches = sum(1 for a, b in zip(s1, s2) if a == b)
    return matches / len(s1)


def _sliding_match(query: str, target: str, window: int, threshold: float) -> list:
    """滑窗比对，找到所有高相似度区段"""
    matches = []
    if len(target) < window:
        return matches
    for i in range(len(query) - window + 1):
        segment = query[i:i + window]
        sim = _hamming_similarity(segment, target[:window])
        if sim >= threshold:
            matches.append({
                "position": i,
                "query_segment": segment,
                "identity": round(sim, 4),
            })
    return matches


def screen_sequence(seq: str) -> dict:
    """
    筛查序列是否包含管制病原体特征序列。

    Args:
        seq: DNA 序列

    Returns:
        dict: 筛查结果

    Raises:
        BiosafetyViolation: 当检测到 Risk Group 4 病原体时
    """
    seq = seq.upper().strip()

    with open(os.path.join(DATA_DIR, "pathogen_signatures.json"), "r") as f:
        pathogen_db = json.load(f)

    config = pathogen_db["screening_config"]
    min_len = config["min_match_length"]
    threshold = config["identity_threshold"]

    alerts = []
    critical = False

    for pathogen in pathogen_db["regulated_pathogens"]:
        name = pathogen["name"]
        risk_group = pathogen["risk_group"]

        for marker in pathogen["markers"]:
            sig = marker["signature"].upper()
            gene = marker["gene"]

            # 正向比对
            matches = _sliding_match(seq, sig, min(len(sig), min_len), threshold)

            # 反向互补比对
            rc_sig = _reverse_complement(sig)
            rc_matches = _sliding_match(seq, rc_sig, min(len(rc_sig), min_len), threshold)

            all_matches = matches + [
                {**m, "strand": "reverse_complement"} for m in rc_matches
            ]

            if all_matches:
                best = max(all_matches, key=lambda m: m["identity"])
                alert = {
                    "pathogen": name,
                    "risk_group": risk_group,
                    "gene": gene,
                    "description": marker["description"],
                    "best_match_identity": best["identity"],
                    "match_position": best["position"],
                    "alert_level": config["alert_level"].get(
                        f"risk_group_{risk_group}", "UNKNOWN"
                    ),
                }
                alerts.append(alert)

                if risk_group >= 4:
                    critical = True

    result = {
        "sequence_length": len(seq),
        "alerts_count": len(alerts),
        "alerts": alerts,
        "safe": len(alerts) == 0,
        "critical": critical,
    }

    if critical:
        # Risk Group 4 → 抛出异常终止
        rg4 = [a for a in alerts if a["risk_group"] >= 4]
        if rg4:
            a = rg4[0]
            raise BiosafetyViolation(
                a["pathogen"], a["risk_group"], a["gene"], a["best_match_identity"]
            )

    return result


def check_virulence_factors(seq: str) -> dict:
    """
    扫描常见毒力因子模式。

    Args:
        seq: DNA 序列
    """
    seq = seq.upper().strip()

    # 简化版毒力因子特征
    virulence_patterns = {
        "Type III Secretion Signal": {
            "pattern": "ATGAACAATAA",  # 简化的 T3SS signal
            "description": "Potential Type III secretion signal peptide",
        },
        "Toxin N-terminal": {
            "pattern": "ATGGATCC",  # 简化
            "description": "Generic toxin-like N-terminal motif",
        },
    }

    findings = []
    for name, info in virulence_patterns.items():
        pattern = info["pattern"]
        pos = 0
        while True:
            idx = seq.find(pattern, pos)
            if idx == -1:
                break
            findings.append({
                "factor": name,
                "position": idx,
                "description": info["description"],
            })
            pos = idx + 1

    return {
        "virulence_factors_found": len(findings),
        "findings": findings[:20],
        "note": "This is a simplified screening. For comprehensive analysis, use VFDB or PHI-base.",
    }


def _reverse_complement(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp.get(b, "N") for b in reversed(seq))
