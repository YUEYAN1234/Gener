"""
结构分析器 - Structure Analyzer
G-quadruplex, Z-DNA, i-Motif, H-DNA 检测
弯曲刚性、结构风险图谱
"""

import re
from typing import Optional


def find_g_quadruplex(seq: str) -> dict:
    """
    扫描 G-quadruplex (G4) 基序。
    模式: G{3,}N{1,7}G{3,}N{1,7}G{3,}N{1,7}G{3,}
    
    Args:
        seq: DNA 序列
    Returns:
        dict: G4 位点列表及风险评估
    """
    seq = seq.upper().strip()
    # 标准 G4 正则
    g4_pattern = re.compile(r'(G{3,6})(\w{1,7})(G{3,6})(\w{1,7})(G{3,6})(\w{1,7})(G{3,6})')
    
    g4_sites = []
    for m in g4_pattern.finditer(seq):
        g_runs = [m.group(1), m.group(3), m.group(5), m.group(7)]
        loops = [m.group(2), m.group(4), m.group(6)]
        min_g = min(len(g) for g in g_runs)
        
        # G4 稳定性评分
        stability = "high" if min_g >= 4 else "moderate" if min_g >= 3 else "low"
        loop_penalty = sum(1 for l in loops if len(l) > 4)
        if loop_penalty >= 2:
            stability = "low"
            
        g4_sites.append({
            "position": m.start(),
            "end": m.end(),
            "sequence": m.group(),
            "length": m.end() - m.start(),
            "g_runs": g_runs,
            "loops": loops,
            "min_g_tract": min_g,
            "stability": stability,
            "synthesis_risk": "HIGH" if stability == "high" else "MODERATE",
        })
    
    # 互补链也检查 (C-rich → G4 on complementary)
    c4_pattern = re.compile(r'(C{3,6})(\w{1,7})(C{3,6})(\w{1,7})(C{3,6})(\w{1,7})(C{3,6})')
    for m in c4_pattern.finditer(seq):
        g4_sites.append({
            "position": m.start(),
            "end": m.end(),
            "sequence": m.group(),
            "length": m.end() - m.start(),
            "strand": "complementary",
            "stability": "moderate",
            "synthesis_risk": "MODERATE",
            "note": "G4 on complementary strand",
        })
    
    return {
        "g4_count": len(g4_sites),
        "g4_sites": g4_sites,
        "overall_risk": "HIGH" if any(s.get("stability") == "high" for s in g4_sites) else
                       "MODERATE" if g4_sites else "LOW",
    }


def find_z_dna(seq: str, min_length: int = 10) -> dict:
    """
    检测 Z-DNA 倾向区段（交替嘌呤-嘧啶序列）。
    
    Args:
        seq: DNA 序列
        min_length: 最小连续交替长度
    """
    seq = seq.upper().strip()
    purine = set("AG")
    pyrimidine = set("CT")
    
    z_regions = []
    i = 0
    while i < len(seq) - 1:
        # 检测交替 purine-pyrimidine
        start = i
        j = i
        while j < len(seq) - 1:
            curr_pur = seq[j] in purine
            next_pur = seq[j + 1] in purine
            if curr_pur != next_pur:  # 交替
                j += 1
            else:
                break
        
        length = j - start + 1
        if length >= min_length:
            segment = seq[start:start + length]
            # GC-rich 的交替序列更容易形成 Z-DNA
            gc_frac = sum(1 for b in segment if b in "GC") / length
            propensity = "high" if gc_frac > 0.6 else "moderate" if gc_frac > 0.4 else "low"
            
            z_regions.append({
                "position": start,
                "end": start + length,
                "sequence": segment,
                "length": length,
                "gc_fraction": round(gc_frac, 3),
                "z_propensity": propensity,
            })
        i = max(j, i + 1)
    
    return {
        "z_dna_count": len(z_regions),
        "z_dna_regions": z_regions,
        "overall_risk": "HIGH" if any(r["z_propensity"] == "high" for r in z_regions) else
                       "MODERATE" if z_regions else "LOW",
    }


def find_i_motif(seq: str) -> dict:
    """
    检测 i-Motif（C-rich 四链体结构）。
    模式: C{3,}N{1,7}C{3,}N{1,7}C{3,}N{1,7}C{3,}
    """
    seq = seq.upper().strip()
    pattern = re.compile(r'(C{3,6})(\w{1,7})(C{3,6})(\w{1,7})(C{3,6})(\w{1,7})(C{3,6})')
    
    sites = []
    for m in pattern.finditer(seq):
        c_runs = [m.group(1), m.group(3), m.group(5), m.group(7)]
        min_c = min(len(c) for c in c_runs)
        sites.append({
            "position": m.start(),
            "end": m.end(),
            "sequence": m.group(),
            "min_c_tract": min_c,
            "stability": "high" if min_c >= 4 else "moderate",
            "note": "Forms at acidic pH (<6.5), may affect in vitro reactions",
        })
    
    return {"i_motif_count": len(sites), "i_motif_sites": sites}


def find_h_dna(seq: str, min_mirror: int = 10) -> dict:
    """
    检测 H-DNA（三股螺旋，mirror repeat）。
    寻找 mirror repeat 序列（正向回文）。
    """
    seq = seq.upper().strip()
    n = len(seq)
    h_sites = []
    
    for length in range(min_mirror, min(30, n // 2)):
        for i in range(n - 2 * length):
            seg1 = seq[i:i + length]
            # mirror repeat: same sequence downstream
            for gap in range(0, min(50, n - i - 2 * length)):
                j = i + length + gap
                if j + length > n:
                    break
                seg2 = seq[j:j + length]
                if seg1 == seg2:
                    # 检查是否富含 purine 或 pyrimidine
                    pur_frac = sum(1 for b in seg1 if b in "AG") / length
                    h_sites.append({
                        "position": i,
                        "mirror_position": j,
                        "sequence": seg1,
                        "gap": gap,
                        "length": length,
                        "purine_fraction": round(pur_frac, 3),
                        "h_dna_propensity": "high" if pur_frac > 0.8 or pur_frac < 0.2 else "moderate",
                    })
            if len(h_sites) > 100:
                break
        if len(h_sites) > 100:
            break
    
    return {
        "h_dna_count": len(h_sites),
        "h_dna_sites": h_sites[:50],
        "overall_risk": "HIGH" if any(s["h_dna_propensity"] == "high" for s in h_sites) else
                       "LOW" if not h_sites else "MODERATE",
    }


def calc_persistence_length(seq: str) -> dict:
    """
    估算 DNA 弯曲刚性（Persistence Length）。
    基于二核苷酸 step 参数。典型值：~50nm (150bp) for B-DNA。
    
    Args:
        seq: DNA 序列
    """
    seq = seq.upper().strip()
    # 二核苷酸弯曲角度（度），基于 Bolshoy et al. 1991
    bend_angles = {
        "AA": 0.0, "AT": 2.0, "AC": 4.7, "AG": 4.0,
        "TA": 6.0, "TT": 0.0, "TC": 4.0, "TG": 7.2,
        "CA": 7.2, "CT": 4.0, "CC": 3.7, "CG": 5.9,
        "GA": 4.7, "GT": 4.7, "GC": 3.0, "GG": 3.7,
    }
    
    angles = []
    for i in range(len(seq) - 1):
        di = seq[i:i+2]
        angles.append(bend_angles.get(di, 4.0))
    
    if not angles:
        return {"error": "Sequence too short"}
    
    import statistics
    avg_angle = statistics.mean(angles)
    # 高弯曲区域
    flexible_regions = []
    window = 20
    for i in range(0, len(angles) - window):
        w_avg = statistics.mean(angles[i:i+window])
        if w_avg > 5.0:
            flexible_regions.append({"position": i, "avg_bend_angle": round(w_avg, 2)})
    
    return {
        "average_bend_angle_deg": round(avg_angle, 2),
        "flexibility": "high" if avg_angle > 4.5 else "normal" if avg_angle > 3.0 else "rigid",
        "estimated_persistence_length_nm": round(50.0 * (3.5 / max(avg_angle, 0.1)), 1),
        "flexible_regions": flexible_regions[:30],
        "total_flexible_regions": len(flexible_regions),
    }


def structural_risk_map(seq: str) -> dict:
    """
    全序列结构风险热图数据。综合 G4、Z-DNA、i-Motif、H-DNA 分析。
    
    Args:
        seq: DNA 序列
    """
    g4 = find_g_quadruplex(seq)
    zdna = find_z_dna(seq)
    imotif = find_i_motif(seq)
    hdna = find_h_dna(seq)
    persist = calc_persistence_length(seq)
    
    # 汇总风险区
    risk_regions = []
    for site in g4.get("g4_sites", []):
        risk_regions.append({
            "type": "G-quadruplex", "position": site["position"],
            "end": site["end"], "risk": site.get("synthesis_risk", "MODERATE"),
        })
    for site in zdna.get("z_dna_regions", []):
        risk_regions.append({
            "type": "Z-DNA", "position": site["position"],
            "end": site["end"], "risk": "HIGH" if site["z_propensity"] == "high" else "MODERATE",
        })
    for site in imotif.get("i_motif_sites", []):
        risk_regions.append({
            "type": "i-Motif", "position": site["position"],
            "end": site["end"], "risk": "MODERATE",
        })
    for site in hdna.get("h_dna_sites", []):
        risk_regions.append({
            "type": "H-DNA", "position": site["position"],
            "end": site.get("mirror_position", site["position"]) + site["length"],
            "risk": "HIGH" if site["h_dna_propensity"] == "high" else "MODERATE",
        })
    
    risk_regions.sort(key=lambda r: r["position"])
    
    overall_risks = [g4["overall_risk"], zdna["overall_risk"]]
    overall = "HIGH" if "HIGH" in overall_risks else "MODERATE" if "MODERATE" in overall_risks else "LOW"
    
    return {
        "sequence_length": len(seq),
        "overall_structural_risk": overall,
        "risk_regions": risk_regions,
        "summary": {
            "g_quadruplex": g4["g4_count"],
            "z_dna": zdna["z_dna_count"],
            "i_motif": imotif["i_motif_count"],
            "h_dna": hdna["h_dna_count"],
        },
        "persistence_length": persist,
    }
