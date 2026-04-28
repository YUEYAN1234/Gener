"""
工具注册表 - Tool Registry
自动扫描 tools/ 下所有工具函数，生成 DeepSeek function calling 的 tools 参数
"""

import inspect
import json
from typing import Callable, Any

# 工具注册存储
_REGISTRY: dict[str, dict] = {}


def tool(
    name: str,
    description: str,
    parameters: dict,
):
    """
    工具注册装饰器。

    Args:
        name: 工具名称（function calling 用）
        description: 工具描述
        parameters: JSON Schema 格式的参数定义
    """
    def decorator(func: Callable) -> Callable:
        _REGISTRY[name] = {
            "function": func,
            "name": name,
            "description": description,
            "parameters": parameters,
        }
        func._tool_name = name
        return func
    return decorator


def get_all_tools() -> list[dict]:
    """返回所有注册工具的 function calling 格式定义"""
    tools = []
    for name, info in _REGISTRY.items():
        tools.append({
            "type": "function",
            "function": {
                "name": info["name"],
                "description": info["description"],
                "parameters": info["parameters"],
            },
        })
    return tools


def execute_tool(name: str, arguments: dict) -> Any:
    """执行指定工具"""
    if name not in _REGISTRY:
        return {"error": f"Unknown tool: {name}. Available: {list(_REGISTRY.keys())}"}
    
    func = _REGISTRY[name]["function"]
    try:
        result = func(**arguments)
        return result
    except Exception as e:
        return {"error": f"Tool execution error: {type(e).__name__}: {str(e)}"}


def get_tool_names() -> list[str]:
    """返回所有注册工具的名称列表"""
    return list(_REGISTRY.keys())


# ============================================================
# 注册所有工具
# ============================================================

from tools.sequence_validator import validate_sequence, calc_gc_content, find_restriction_sites, basic_stats
from tools.thermodynamics import calc_nn_thermodynamics, calc_tm, predict_hairpins, calc_overlap_compatibility
from tools.structure_analyzer import (
    find_g_quadruplex, find_z_dna, find_i_motif, find_h_dna,
    calc_persistence_length, structural_risk_map,
)
from tools.complexity_scorer import score_complexity, identify_difficult_regions
from tools.loxp_designer import design_loxp_sites, check_orthogonality, simulate_scramble
from tools.telomere_centromere import (
    design_telomere, design_centromere, assess_stability, optimize_for_synthesis,
)
from tools.fragment_assembler import (
    split_gibson, split_golden_gate, optimize_overlaps,
    generate_assembly_report, suggest_pcr_protocol,
)
from tools.biosafety_checker import screen_sequence, check_virulence_factors
from tools.format_converter import to_fasta, to_genbank, export_fragment_list


# --- 序列验证工具 ---

@tool(
    name="validate_sequence",
    description="验证 DNA 序列的合法性。检查序列是否仅包含合法 IUPAC 碱基字符。",
    parameters={
        "type": "object",
        "properties": {
            "seq": {"type": "string", "description": "DNA 序列字符串"},
        },
        "required": ["seq"],
    },
)
def _validate_sequence(seq: str):
    return validate_sequence(seq)


@tool(
    name="calc_gc_content",
    description="计算 DNA 序列的 GC 含量（全局 & 滑窗），识别极端 GC 区段。",
    parameters={
        "type": "object",
        "properties": {
            "seq": {"type": "string", "description": "DNA 序列"},
            "window": {"type": "integer", "description": "滑窗大小 (bp)，0 表示仅全局", "default": 100},
        },
        "required": ["seq"],
    },
)
def _calc_gc_content(seq: str, window: int = 100):
    return calc_gc_content(seq, window)


@tool(
    name="find_restriction_sites",
    description="扫描 DNA 序列中的限制性酶切位点。可指定特定酶列表或扫描全部 40+ 种酶。",
    parameters={
        "type": "object",
        "properties": {
            "seq": {"type": "string", "description": "DNA 序列"},
            "enzymes": {"type": "array", "items": {"type": "string"}, "description": "指定酶列表，如 ['EcoRI','BamHI']。不提供则扫描全部"},
        },
        "required": ["seq"],
    },
)
def _find_restriction_sites(seq: str, enzymes=None):
    return find_restriction_sites(seq, enzymes)


@tool(
    name="basic_stats",
    description="计算 DNA 序列的基本统计信息：长度、碱基组成、AT/GC skew、简单重复序列。",
    parameters={
        "type": "object",
        "properties": {
            "seq": {"type": "string", "description": "DNA 序列"},
        },
        "required": ["seq"],
    },
)
def _basic_stats(seq: str):
    return basic_stats(seq)


# --- 热力学工具 ---

@tool(
    name="calc_nn_thermodynamics",
    description="使用 SantaLucia 1998 最近邻模型计算 DNA 双链热力学参数：ΔH, ΔS, ΔG(37°C), Tm。包含盐浓度和 Mg2+ 校正。公式：ΔG_total = ΔG_initiation + Σ ΔG_stacking",
    parameters={
        "type": "object",
        "properties": {
            "seq": {"type": "string", "description": "DNA 序列（5'→3'）"},
            "na": {"type": "number", "description": "[Na+] 浓度 (mM)", "default": 50},
            "mg": {"type": "number", "description": "[Mg2+] 浓度 (mM)", "default": 1.5},
            "dnac": {"type": "number", "description": "DNA 链总浓度 (nM)", "default": 250},
        },
        "required": ["seq"],
    },
)
def _calc_nn_thermodynamics(seq, na=50.0, mg=1.5, dnac=250.0):
    return calc_nn_thermodynamics(seq, na, mg, dnac)


@tool(
    name="calc_tm",
    description="快速计算 DNA 序列的 Tm 值。短序列(<60bp)用最近邻模型，长序列用 Marmur-Doty 经验公式。默认条件：[Na+]=50mM, [Mg2+]=1.5mM。",
    parameters={
        "type": "object",
        "properties": {
            "seq": {"type": "string", "description": "DNA 序列"},
            "na": {"type": "number", "description": "[Na+] (mM)", "default": 50},
            "mg": {"type": "number", "description": "[Mg2+] (mM)", "default": 1.5},
        },
        "required": ["seq"],
    },
)
def _calc_tm(seq, na=50.0, mg=1.5):
    return calc_tm(seq, na, mg)


@tool(
    name="predict_hairpins",
    description="预测 DNA 序列中的发夹结构 (Hairpin)。识别可能阻碍聚合酶滑动的二级结构。高风险：ΔG < -3 kcal/mol。",
    parameters={
        "type": "object",
        "properties": {
            "seq": {"type": "string", "description": "DNA 序列"},
            "min_stem": {"type": "integer", "description": "最小 stem 长度", "default": 4},
        },
        "required": ["seq"],
    },
)
def _predict_hairpins(seq, min_stem=4):
    return predict_hairpins(seq, min_stem)


@tool(
    name="calc_overlap_compatibility",
    description="检查多个 DNA 片段的 Overlap 区 Tm 均衡性。确保所有 overlap 的 Tm 差值 < 5°C，避免非特异性组装。",
    parameters={
        "type": "object",
        "properties": {
            "fragments": {
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {
                        "name": {"type": "string"},
                        "overlap_3prime": {"type": "string"},
                    },
                },
                "description": "片段列表",
            },
            "target_tm": {"type": "number", "description": "目标 Tm (°C)", "default": 60},
        },
        "required": ["fragments"],
    },
)
def _calc_overlap_compatibility(fragments, target_tm=60.0):
    return calc_overlap_compatibility(fragments, target_tm)


# --- 结构分析工具 ---

@tool(
    name="find_g_quadruplex",
    description="扫描 DNA 序列中的 G-quadruplex (G4) 结构基序。G4 会严重阻碍 PCR 扩增和测序。模式: G3+N1-7G3+N1-7G3+N1-7G3+。",
    parameters={
        "type": "object",
        "properties": {"seq": {"type": "string", "description": "DNA 序列"}},
        "required": ["seq"],
    },
)
def _find_g_quadruplex(seq):
    return find_g_quadruplex(seq)


@tool(
    name="find_z_dna",
    description="检测 Z-DNA 倾向区段。交替嘌呤-嘧啶序列在负超螺旋应力下可形成左旋 Z-DNA。",
    parameters={
        "type": "object",
        "properties": {
            "seq": {"type": "string", "description": "DNA 序列"},
            "min_length": {"type": "integer", "description": "最小连续交替长度", "default": 10},
        },
        "required": ["seq"],
    },
)
def _find_z_dna(seq, min_length=10):
    return find_z_dna(seq, min_length)


@tool(
    name="find_i_motif",
    description="检测 i-Motif 结构（C-rich 四链体）。在酸性 pH 下形成，可能影响体外反应。",
    parameters={
        "type": "object",
        "properties": {"seq": {"type": "string", "description": "DNA 序列"}},
        "required": ["seq"],
    },
)
def _find_i_motif(seq):
    return find_i_motif(seq)


@tool(
    name="find_h_dna",
    description="检测 H-DNA 三股螺旋结构（mirror repeat）。富含嘌呤或嘧啶的 mirror repeat 可形成三链体。",
    parameters={
        "type": "object",
        "properties": {"seq": {"type": "string", "description": "DNA 序列"}},
        "required": ["seq"],
    },
)
def _find_h_dna(seq):
    return find_h_dna(seq)


@tool(
    name="structural_risk_map",
    description="生成全序列结构风险热图。综合 G4、Z-DNA、i-Motif、H-DNA 分析 + 弯曲刚性，给出整体风险评估。",
    parameters={
        "type": "object",
        "properties": {"seq": {"type": "string", "description": "DNA 序列"}},
        "required": ["seq"],
    },
)
def _structural_risk_map(seq):
    return structural_risk_map(seq)


# --- 难度评分 ---

@tool(
    name="score_complexity",
    description="综合合成难度评分 (0-100)。考虑 GC 含量、重复序列、G4/Z-DNA 密度、同源区段等因子。",
    parameters={
        "type": "object",
        "properties": {"seq": {"type": "string", "description": "DNA 序列"}},
        "required": ["seq"],
    },
)
def _score_complexity(seq):
    return score_complexity(seq)


@tool(
    name="identify_difficult_regions",
    description="标记序列中的具体合成困难区段及原因（高GC、G4、Z-DNA 等）。",
    parameters={
        "type": "object",
        "properties": {"seq": {"type": "string", "description": "DNA 序列"}},
        "required": ["seq"],
    },
)
def _identify_difficult_regions(seq):
    return identify_difficult_regions(seq)


# --- loxP 设计 ---

@tool(
    name="design_loxp_sites",
    description="在 DNA 序列中自动布置 loxP 位点用于 SCRaMbLE 进化工程。策略性避开必需基因，在非必需基因下游或基因间区 (IGS) 布点。",
    parameters={
        "type": "object",
        "properties": {
            "seq": {"type": "string", "description": "DNA 序列"},
            "interval": {"type": "integer", "description": "loxP 间隔 (bp)", "default": 5000},
            "variant": {"type": "string", "description": "loxP 变体: loxP, lox511, lox2272", "default": "loxP"},
            "avoid_regions": {
                "type": "array",
                "items": {"type": "object", "properties": {
                    "start": {"type": "integer"}, "end": {"type": "integer"}, "name": {"type": "string"}
                }},
                "description": "需要避开的区域（如必需基因）",
            },
        },
        "required": ["seq"],
    },
)
def _design_loxp_sites(seq, interval=5000, variant="loxP", avoid_regions=None):
    return design_loxp_sites(seq, interval, variant, avoid_regions)


@tool(
    name="check_orthogonality",
    description="验证 loxP 变体间的正交性。确保不同 spacer 的 lox 位点不会交叉重组。",
    parameters={
        "type": "object",
        "properties": {
            "variants": {"type": "array", "items": {"type": "string"}, "description": "要检查的变体列表"},
        },
        "required": [],
    },
)
def _check_orthogonality(variants=None):
    return check_orthogonality(variants)


@tool(
    name="simulate_scramble",
    description="模拟 Cre 酶介导的 SCRaMbLE 拓扑演变。预测多 loxP 系统在重组酶作用下的删除和反转事件。",
    parameters={
        "type": "object",
        "properties": {
            "seq_length": {"type": "integer", "description": "序列总长度"},
            "loxp_positions": {"type": "array", "items": {"type": "integer"}, "description": "loxP 位点位置列表"},
            "n_events": {"type": "integer", "description": "模拟事件数", "default": 100},
        },
        "required": ["seq_length", "loxp_positions"],
    },
)
def _simulate_scramble(seq_length, loxp_positions, n_events=100):
    return simulate_scramble(seq_length, loxp_positions, n_events)


# --- 端粒/着丝粒 ---

@tool(
    name="design_telomere",
    description="生成端粒重复序列。支持多物种（yeast, human, arabidopsis 等），评估 G4 合成风险。",
    parameters={
        "type": "object",
        "properties": {
            "organism": {"type": "string", "description": "目标物种", "default": "yeast"},
            "length": {"type": "integer", "description": "端粒长度 (bp)", "default": 2000},
            "add_subtelomeric": {"type": "boolean", "description": "是否添加亚端粒缓冲区", "default": True},
        },
        "required": [],
    },
)
def _design_telomere(organism="yeast", length=2000, add_subtelomeric=True):
    return design_telomere(organism, length, add_subtelomeric)


@tool(
    name="design_centromere",
    description="设计着丝粒序列。当前支持 S. cerevisiae 点着丝粒 (CDE I/II/III)。",
    parameters={
        "type": "object",
        "properties": {
            "organism": {"type": "string", "description": "目标物种", "default": "yeast"},
            "chromosome_id": {"type": "string", "description": "染色体编号", "default": "synIII"},
        },
        "required": [],
    },
)
def _design_centromere(organism="yeast", chromosome_id="synIII"):
    return design_centromere(organism, chromosome_id)


@tool(
    name="assess_stability",
    description="评估功能元件（端粒、着丝粒等）在宿主中的复制稳定性。",
    parameters={
        "type": "object",
        "properties": {
            "element_seq": {"type": "string", "description": "元件 DNA 序列"},
            "host": {"type": "string", "description": "宿主物种", "default": "yeast"},
        },
        "required": ["element_seq"],
    },
)
def _assess_stability(element_seq, host="yeast"):
    return assess_stability(element_seq, host)


# --- 片段组装 ---

@tool(
    name="split_gibson",
    description="使用 Gibson Assembly 策略自动拆分大片段 DNA。计算 Overlap 序列及退火温度，优化切割位点。",
    parameters={
        "type": "object",
        "properties": {
            "seq": {"type": "string", "description": "DNA 序列"},
            "target_size": {"type": "integer", "description": "目标片段大小 (bp)", "default": 5000},
            "overlap": {"type": "integer", "description": "Overlap 长度 (bp)", "default": 40},
        },
        "required": ["seq"],
    },
)
def _split_gibson(seq, target_size=5000, overlap=40):
    return split_gibson(seq, target_size, overlap)


@tool(
    name="split_golden_gate",
    description="使用 Golden Gate Assembly 策略拆分序列。利用 Type IIS 限制酶生成唯一的 4bp overhang。",
    parameters={
        "type": "object",
        "properties": {
            "seq": {"type": "string", "description": "DNA 序列"},
            "target_size": {"type": "integer", "description": "目标片段大小 (bp)", "default": 5000},
            "enzyme": {"type": "string", "description": "Type IIS 酶", "default": "BsaI"},
        },
        "required": ["seq"],
    },
)
def _split_golden_gate(seq, target_size=5000, enzyme="BsaI"):
    return split_golden_gate(seq, target_size, enzyme)


@tool(
    name="generate_assembly_report",
    description="生成完整的 DNA 组装方案报告，包含片段信息、实验协议建议。",
    parameters={
        "type": "object",
        "properties": {
            "fragments": {"type": "array", "items": {"type": "object"}, "description": "片段列表"},
            "method": {"type": "string", "description": "组装方法: Gibson / Golden Gate", "default": "Gibson"},
        },
        "required": ["fragments"],
    },
)
def _generate_assembly_report(fragments, method="Gibson"):
    return generate_assembly_report(fragments, method)


@tool(
    name="suggest_pcr_protocol",
    description="为特定 DNA 片段建议 PCR 扩增条件。针对高 GC、长片段等困难模板自动调整参数和添加剂。",
    parameters={
        "type": "object",
        "properties": {
            "fragment": {"type": "object", "description": "片段信息 {name, length, gc_content}",
                        "properties": {
                            "name": {"type": "string"},
                            "length": {"type": "integer"},
                            "gc_content": {"type": "number"},
                        }},
        },
        "required": ["fragment"],
    },
)
def _suggest_pcr_protocol(fragment):
    return suggest_pcr_protocol(fragment)


# --- 生物安全 ---

@tool(
    name="screen_sequence",
    description="🚨 生物安全筛查：检查 DNA 序列是否包含管制病原体（埃博拉、天花、炭疽等）的特征序列。Risk Group 4 触发立即终止。",
    parameters={
        "type": "object",
        "properties": {"seq": {"type": "string", "description": "DNA 序列"}},
        "required": ["seq"],
    },
)
def _screen_sequence(seq):
    return screen_sequence(seq)


# --- 格式转换 ---

@tool(
    name="to_fasta",
    description="将 DNA 序列转换为 FASTA 格式文本。",
    parameters={
        "type": "object",
        "properties": {
            "seq": {"type": "string", "description": "DNA 序列"},
            "name": {"type": "string", "description": "序列名称", "default": "synthetic_construct"},
            "description": {"type": "string", "description": "序列描述", "default": ""},
        },
        "required": ["seq"],
    },
)
def _to_fasta(seq, name="synthetic_construct", description=""):
    return to_fasta(seq, name, description)


@tool(
    name="to_genbank",
    description="将 DNA 序列转换为 GenBank 格式文本，包含 Feature Table。",
    parameters={
        "type": "object",
        "properties": {
            "seq": {"type": "string", "description": "DNA 序列"},
            "name": {"type": "string", "description": "序列名称", "default": "synthetic_construct"},
            "organism": {"type": "string", "description": "来源物种", "default": "synthetic construct"},
            "features": {
                "type": "array",
                "items": {"type": "object"},
                "description": "特征列表 [{type, start, end, qualifiers: {}}]",
            },
        },
        "required": ["seq"],
    },
)
def _to_genbank(seq, name="synthetic_construct", organism="synthetic construct", features=None):
    return to_genbank(seq, name, organism, features)


@tool(
    name="export_fragment_list",
    description="导出组装片段清单。支持 text / csv / json 格式。",
    parameters={
        "type": "object",
        "properties": {
            "fragments": {"type": "array", "items": {"type": "object"}, "description": "片段列表"},
            "format": {"type": "string", "description": "输出格式: text / csv / json", "default": "text"},
        },
        "required": ["fragments"],
    },
)
def _export_fragment_list(fragments, format="text"):
    return export_fragment_list(fragments, format)


# --- NCBI 数据库检索 ---

from tools.ncbi_search import (
    search_nucleotide, fetch_sequence, search_gene,
    fetch_gene_sequence, blast_short,
)


@tool(
    name="search_nucleotide",
    description="在 NCBI Nucleotide 数据库中搜索 DNA/RNA 序列。可按基因名、功能描述、Accession 等关键词检索，支持限定物种。返回 accession、描述、长度等摘要信息。",
    parameters={
        "type": "object",
        "properties": {
            "query": {"type": "string", "description": "搜索关键词（基因名、功能描述、accession 等）"},
            "max_results": {"type": "integer", "description": "最大返回条数 (1-20)", "default": 10},
            "organism": {"type": "string", "description": "限定物种，如 'Saccharomyces cerevisiae'"},
        },
        "required": ["query"],
    },
)
def _search_nucleotide(query, max_results=10, organism=None):
    return search_nucleotide(query, max_results, organism)


@tool(
    name="fetch_sequence",
    description="通过 NCBI Accession Number 获取完整的核酸序列及注释。支持 FASTA 和 GenBank 格式，可指定截取范围。用于获取参考序列、已有基因组片段等。",
    parameters={
        "type": "object",
        "properties": {
            "accession": {"type": "string", "description": "NCBI Accession Number（如 NC_001133.9, NM_001301302）"},
            "seq_start": {"type": "integer", "description": "截取起始位置 (1-based，可选)"},
            "seq_end": {"type": "integer", "description": "截取终止位置 (1-based，可选)"},
            "format": {"type": "string", "description": "返回格式: fasta 或 genbank", "default": "fasta"},
        },
        "required": ["accession"],
    },
)
def _fetch_sequence(accession, seq_start=None, seq_end=None, format="fasta"):
    return fetch_sequence(accession, seq_start, seq_end, format)


@tool(
    name="search_gene",
    description="在 NCBI Gene 数据库中搜索基因信息。返回基因名、描述、物种、染色体位置、功能摘要等。适用于查找目标基因的详细注释和定位信息。",
    parameters={
        "type": "object",
        "properties": {
            "query": {"type": "string", "description": "基因名或关键词"},
            "organism": {"type": "string", "description": "限定物种"},
            "max_results": {"type": "integer", "description": "最大返回条数", "default": 10},
        },
        "required": ["query"],
    },
)
def _search_gene(query, organism=None, max_results=10):
    return search_gene(query, organism, max_results)


@tool(
    name="fetch_gene_sequence",
    description="通过 NCBI Gene ID 获取基因关联的核酸序列列表（mRNA、基因组序列等 Accession）。可配合 fetch_sequence 使用获取完整序列。",
    parameters={
        "type": "object",
        "properties": {
            "gene_id": {"type": "string", "description": "NCBI Gene ID（纯数字）"},
        },
        "required": ["gene_id"],
    },
)
def _fetch_gene_sequence(gene_id):
    return fetch_gene_sequence(gene_id)


@tool(
    name="blast_short",
    description="对短片段 DNA 序列（<1000bp）执行 NCBI BLAST 在线比对。用于检查序列同源性、确认序列来源。注意：耗时较长(30-120秒)，仅建议短片段使用。",
    parameters={
        "type": "object",
        "properties": {
            "sequence": {"type": "string", "description": "查询 DNA 序列 (<1000bp)"},
            "database": {"type": "string", "description": "数据库: nt(核酸) 或 refseq_rna", "default": "nt"},
            "max_hits": {"type": "integer", "description": "最大命中数", "default": 5},
        },
        "required": ["sequence"],
    },
)
def _blast_short(sequence, database="nt", max_hits=5):
    return blast_short(sequence, database, max_hits)

