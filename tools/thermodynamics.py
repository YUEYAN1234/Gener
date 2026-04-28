"""
热力学引擎 - Thermodynamics Engine
最近邻模型计算 ΔG, ΔH, ΔS, Tm
发夹结构预测、Overlap 兼容性检查
"""

import json
import math
import os
import statistics
from typing import Optional

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

with open(os.path.join(DATA_DIR, "nn_parameters.json"), "r") as f:
    NN_PARAMS = json.load(f)

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
R = 1.987  # cal/(mol·K)


def _complement(seq):
    return "".join(COMPLEMENT.get(b, "N") for b in seq)

def _reverse_complement(seq):
    return _complement(seq)[::-1]

def _is_self_complementary(seq):
    return seq == _reverse_complement(seq)


def calc_nn_thermodynamics(seq: str, na: float = 50.0, mg: float = 1.5, dnac: float = 250.0) -> dict:
    """
    最近邻模型计算 DNA 双链热力学参数。
    ΔG_total = ΔG_initiation + Σ ΔG_stacking + ΔG_penalty

    Args:
        seq: DNA 序列（5'→3'）
        na: [Na+] mM, mg: [Mg2+] mM, dnac: DNA浓度 nM
    """
    seq = seq.upper().strip()
    if len(seq) < 2:
        return {"error": "Sequence must be >= 2 bases"}

    stacking = NN_PARAMS["stacking"]
    init_params = NN_PARAMS["initiation"]
    dH, dS = 0.0, 0.0

    # Initiation
    for end_base in [seq[0], seq[-1]]:
        key = "init_gc" if end_base in "GC" else "init_at"
        dH += init_params[key]["dH"]
        dS += init_params[key]["dS"]

    # Stacking
    for i in range(len(seq) - 1):
        dt = seq[i:i+2]
        db = _reverse_complement(dt)
        key1 = f"{dt}/{db}"
        key2 = f"{db[::-1]}/{dt[::-1]}"
        if key1 in stacking:
            dH += stacking[key1]["dH"]; dS += stacking[key1]["dS"]
        elif key2 in stacking:
            dH += stacking[key2]["dH"]; dS += stacking[key2]["dS"]
        else:
            dH += -8000; dS += -22.0  # fallback

    if _is_self_complementary(seq):
        dS += NN_PARAMS["symmetry_correction"]["dS"]

    dG_37 = dH - 310.15 * dS
    ct = dnac * 1e-9
    denom = dS + R * math.log(ct if _is_self_complementary(seq) else ct / 4)
    tm_1m = (dH / denom - 273.15) if denom != 0 else 0

    # Salt correction (Owczarzy 2004)
    fgc = sum(1 for b in seq if b in "GC") / len(seq)
    na_m = na / 1000.0
    tm_salt = tm_1m
    if na_m > 0:
        owc = NN_PARAMS["salt_correction"]["owczarzy_2004"]
        inv = 1.0 / (tm_1m + 273.15) + (owc["a"] * fgc + owc["b"]) * math.log(na_m) + owc["c"] * math.log(na_m) ** 2
        tm_salt = 1.0 / inv - 273.15

    # Mg correction (Owczarzy 2008)
    mg_m = mg / 1000.0
    tm_mg = tm_salt
    if mg_m > 0:
        mc = NN_PARAMS["mg_correction"]["coefficients"]
        ln_mg = math.log(mg_m)
        inv = (1.0 / (tm_1m + 273.15) + mc["a"] + mc["b"] * ln_mg
               + fgc * (mc["c"] + mc["d"] * ln_mg)
               + (1.0 / (2 * max(len(seq) - 1, 1))) * (mc["e"] + mc["f"] * ln_mg + mc["g"] * ln_mg ** 2))
        tm_mg = 1.0 / inv - 273.15

    return {
        "sequence_length": len(seq), "dH_cal_mol": round(dH, 1),
        "dS_cal_mol_K": round(dS, 2), "dG_37_cal_mol": round(dG_37, 1),
        "dG_37_kcal_mol": round(dG_37 / 1000, 2),
        "Tm_1M_NaCl_C": round(tm_1m, 1),
        "Tm_salt_corrected_C": round(tm_salt, 1),
        "Tm_Mg_corrected_C": round(tm_mg, 1),
        "conditions": {"Na_mM": na, "Mg_mM": mg, "DNA_conc_nM": dnac},
    }


def calc_tm(seq: str, na: float = 50.0, mg: float = 1.5, dnac: float = 250.0) -> dict:
    """快速 Tm 计算。短序列用 NN，长序列用经验公式。"""
    seq = seq.upper().strip()
    n = len(seq)
    if n < 6:
        return {"error": "Sequence too short"}
    if n <= 60:
        r = calc_nn_thermodynamics(seq, na, mg, dnac)
        return {"Tm_C": r["Tm_Mg_corrected_C"], "method": "Nearest-Neighbor (SantaLucia 1998)", "conditions": r["conditions"]}
    fgc = sum(1 for b in seq if b in "GC") / n
    tm = 81.5 + 16.6 * math.log10(na / 1000.0) + 41.0 * fgc - 600.0 / n
    return {"Tm_C": round(tm, 1), "method": "Marmur-Doty (empirical)", "conditions": {"Na_mM": na, "Mg_mM": mg}}


def predict_hairpins(seq: str, min_stem: int = 4, min_loop: int = 3, max_loop: int = 8) -> dict:
    """预测发夹结构。返回按 ΔG 排序的 hairpin 列表。"""
    seq = seq.upper().strip()
    n = len(seq)
    hairpins = []

    for i in range(n):
        for stem_len in range(min_stem, min(15, (n - i) // 2)):
            for loop_len in range(min_loop, min(max_loop + 1, n - i - stem_len)):
                j = i + stem_len + loop_len
                if j + stem_len > n:
                    break
                stem5 = seq[i:i + stem_len]
                stem3 = seq[j:j + stem_len]
                rc = _reverse_complement(stem3)
                mm = sum(1 for a, b in zip(stem5, rc) if a != b)
                if mm <= max(1, stem_len * 0.15):
                    thermo = calc_nn_thermodynamics(stem5)
                    dg = thermo.get("dG_37_kcal_mol", 0)
                    hairpins.append({
                        "position": i, "stem_5prime": stem5,
                        "loop": seq[i + stem_len:j], "stem_3prime": stem3,
                        "stem_length": stem_len, "loop_length": loop_len,
                        "mismatches": mm, "dG_kcal_mol": dg,
                    })
        if len(hairpins) > 200:
            break

    hairpins.sort(key=lambda h: h["dG_kcal_mol"])
    high_risk = [h for h in hairpins if h["dG_kcal_mol"] < -3.0]
    return {"total_hairpins": len(hairpins), "hairpins": hairpins[:50],
            "high_risk_count": len(high_risk), "high_risk": high_risk[:20]}


def calc_overlap_compatibility(fragments: list, target_tm: float = 60.0, na: float = 50.0, mg: float = 1.5) -> dict:
    """检查多片段 Overlap Tm 均衡性。"""
    overlaps, tms = [], []
    for i, frag in enumerate(fragments):
        ol = frag.get("overlap_3prime", "")
        if ol:
            tm = calc_tm(ol, na, mg).get("Tm_C", 0)
            tms.append(tm)
            next_name = fragments[i + 1].get("name", f"F{i+1}") if i + 1 < len(fragments) else "END"
            overlaps.append({
                "fragment_pair": f"{frag.get('name', f'F{i}')}-{next_name}",
                "overlap_sequence": ol, "overlap_length": len(ol),
                "Tm_C": tm, "deviation_from_target": round(tm - target_tm, 1),
            })
    if not tms:
        return {"error": "No overlaps found"}
    tm_range = max(tms) - min(tms)
    return {
        "overlaps": overlaps,
        "tm_mean_C": round(statistics.mean(tms), 1),
        "tm_std_C": round(statistics.stdev(tms) if len(tms) > 1 else 0, 2),
        "tm_range_C": round(tm_range, 1),
        "compatible": tm_range < 5.0,
        "recommendation": "All overlaps compatible." if tm_range < 5.0
            else f"Tm range ({tm_range:.1f}°C) > 5°C. Adjust overlap lengths.",
    }
