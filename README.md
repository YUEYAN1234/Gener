# 大片段 DNA 进化设计与物化分析智能体 Gener
**Synthetic Biology AI Agent for Large-Fragment DNA Design**

Gener是一个专为合成生物学与生物信息学架构设计的智能体（Agent）系统。系统以 DeepSeek V4-Pro API 作为核心推理引擎，通过 function calling 技术深度整合了 34 个专业的计算生物学与物理化学分析工具，旨在解决大片段 DNA（>50kb）设计、人工染色体构建及基因组进化工程中的复杂问题。

## 核心功能与特性

- **基于 DeepSeek API 的智能推理**：借助先进的大语言模型，Agent 能够理解复杂的自然语言指令，并将其转化为结构化的生物学设计工作流。
- **物理化学预测引擎**：内置基于最近邻模型（SantaLucia 1998）的热力学引擎，支持精确的 $\Delta G$、$\Delta H$、$\Delta S$ 及 $T_m$ 计算，并提供盐离子（$Na^+$）和镁离子（$Mg^{2+}$）浓度校正。
- **高级结构动力学分析**：支持 G-quadruplex (G4)、Z-DNA、i-Motif 及 H-DNA 等非经典拓扑结构的精准定位，并能评估 DNA 弯曲刚性与局部合成风险。
- **数据库检索与序列获取**：集成 NCBI Entrez 接口，支持自动化搜索并提取指定的 Nucleotide (核酸) 和 Gene (基因) 序列及其注释。
- **功能元件与进化工程设计**：
  - 自动设计端粒（Telomeres）和着丝粒（Centromeres）元件。
  - 智能锚定 loxP 位点及其正交变体（lox511, lox2272），支持 Cre 介导的 SCRaMbLE 拓扑演变模拟。
- **大片段组装与拆解**：支持基于 Gibson Assembly 和 Golden Gate 策略的序列自动化拆解，智能优化 Overlap 区域的热力学兼容性。
- **生物安全合规控制**：内置针对受管制病原体特征序列的强制性生物安全筛查机制，触发生物安全红线（如 Risk Group 4 病原体）时立即抛出异常并终止任务。

## 项目架构

项目采用模块化设计，主要包含以下核心组件：

- **`main.py`**: 系统的命令行（CLI）交互入口。
- **`agent/`**: 包含智能体核心引擎（`core.py`）、DeepSeek 客户端封装（`deepseek_client.py`）、系统角色指令（`system_prompt.py`）以及基于反射的工具注册表（`tool_registry.py`）。
- **`tools/`**: 10 个领域的专业生物学分析工具包，暴露 34 个可供 LLM 调用的函数。
- **`data/`**: 包含热力学参数、限制性内切酶数据库、病原体特征序列以及多物种密码子频率表。

## 环境依赖与安装

### 系统要求
- **操作系统**: Windows, Linux, 或 macOS
- **Python 环境**: Python 3.12 或更高版本

### 安装步骤

1. **克隆或下载项目**至本地目录。
2. **创建并激活虚拟环境**：
   ```bash
   python -m venv venv
   # Windows:
   .\venv\Scripts\activate
   # Linux/macOS:
   source venv/bin/activate
   ```
3. **安装依赖包**：
   ```bash
   pip install -r requirements.txt
   ```

## 配置说明

在项目根目录下，修改或创建 `.env` 文件，配置 DeepSeek API 的相关凭证：

```ini
DEEPSEEK_API_KEY=your_deepseek_api_key_here
DEEPSEEK_BASE_URL=https://api.deepseek.com
DEEPSEEK_MODEL=deepseek-chat
```
*注：本项目设计依赖于支持 Function Calling 的模型版本。请确保模型名为 `deepseek-chat`。*

## 运行指南

激活虚拟环境后，通过以下命令启动交互式命令行终端：

```bash
python main.py
```

### 终端内可用命令：
- `/tools` - 列出系统当前加载的所有可用计算工具
- `/reset` - 清除当前对话上下文并重置 Agent 状态
- `/help`  - 显示使用帮助与示例
- `/quit` 或 `/exit` - 安全退出系统

### 交互示例
您可以直接使用自然语言描述复杂的生物学任务，例如：
> "帮我从 NCBI 搜索人类的 p53 基因，提取其序列，然后评估它的 GC 含量和 G4 结构风险。"

> "我想设计一段 80kb 的合成酵母染色体臂，包含一个 2kb 的端粒末端，并且需要在每 5kb 的非必需基因区加入 loxP 位点以便进化。请分析其物理化学稳定性并拆分成 5kb 的合成片段。"

## 工具集概览 (Tools Inventory)

系统共注册了 34 个工具，分类如下：
1. **序列基础处理**: `validate_sequence`, `calc_gc_content`, `find_restriction_sites`, `basic_stats`
2. **数据库检索(NCBI)**: `search_nucleotide`, `fetch_sequence`, `search_gene`, `fetch_gene_sequence`, `blast_short`
3. **热力学分析**: `calc_nn_thermodynamics`, `calc_tm`, `predict_hairpins`, `calc_overlap_compatibility`
4. **拓扑结构分析**: `find_g_quadruplex`, `find_z_dna`, `find_i_motif`, `find_h_dna`, `structural_risk_map`
5. **合成难度评估**: `score_complexity`, `identify_difficult_regions`
6. **进化工程与 loxP**: `design_loxp_sites`, `check_orthogonality`, `simulate_scramble`
7. **染色体元件设计**: `design_telomere`, `design_centromere`, `assess_stability`
8. **片段组装策略**: `split_gibson`, `split_golden_gate`, `generate_assembly_report`, `suggest_pcr_protocol`
9. **安全筛查与导出**: `screen_sequence`, `to_fasta`, `to_genbank`, `export_fragment_list`

## 生物安全声明

本项目严格遵循合成生物学研究的安全伦理规范。系统内置的 `screen_sequence` 工具强制对所有输入与生成的序列进行致病因子检查。任何涉嫌组装受管制的高危病原体（如埃博拉、天花等）的指令将直接触发警报并终止运行过程。

## 许可证

[MIT License] - 仅限学术研究及合规的工业研发用途。
