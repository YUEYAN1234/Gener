"""
格式转换器 - Format Converter
GenBank / FASTA 格式输出，片段清单导出
"""

from datetime import datetime
from typing import Optional


def to_fasta(seq: str, name: str = "synthetic_construct", description: str = "") -> dict:
    """
    生成 FASTA 格式文本。

    Args:
        seq: DNA 序列
        name: 序列名称
        description: 描述
    """
    seq = seq.upper().strip()
    header = f">{name} {description}".strip()

    # 每行 70 字符
    lines = [header]
    for i in range(0, len(seq), 70):
        lines.append(seq[i:i + 70])

    fasta_text = "\n".join(lines)

    return {
        "format": "FASTA",
        "name": name,
        "sequence_length": len(seq),
        "content": fasta_text,
    }


def to_genbank(
    seq: str,
    name: str = "synthetic_construct",
    organism: str = "synthetic construct",
    features: Optional[list] = None,
    molecule_type: str = "DNA",
    topology: str = "linear",
) -> dict:
    """
    生成 GenBank 格式文本。

    Args:
        seq: DNA 序列
        name: 序列名称
        organism: 来源物种
        features: 特征列表 [{type, start, end, qualifiers: {key: value}}]
        molecule_type: DNA / RNA
        topology: linear / circular
    """
    seq = seq.upper().strip()
    n = len(seq)
    features = features or []
    date_str = datetime.now().strftime("%d-%b-%Y").upper()
    locus_name = name[:16].replace(" ", "_")

    lines = []

    # LOCUS line
    lines.append(
        f"LOCUS       {locus_name:<16} {n:>6} bp    {molecule_type:<6} {topology:<8} {date_str}"
    )
    lines.append(f"DEFINITION  {name}.")
    lines.append(f"ACCESSION   .")
    lines.append(f"VERSION     .")
    lines.append(f"SOURCE      {organism}")
    lines.append(f"  ORGANISM  {organism}")
    lines.append(f"            Unclassified.")

    # FEATURES
    lines.append("FEATURES             Location/Qualifiers")
    lines.append(f'     source          1..{n}')
    lines.append(f'                     /organism="{organism}"')
    lines.append(f'                     /mol_type="genomic {molecule_type}"')

    for feat in features:
        feat_type = feat.get("type", "misc_feature")
        start = feat.get("start", 1)
        end = feat.get("end", n)
        strand = feat.get("strand", "+")

        if strand == "-":
            location = f"complement({start}..{end})"
        else:
            location = f"{start}..{end}"

        lines.append(f"     {feat_type:<16}{location}")

        for key, value in feat.get("qualifiers", {}).items():
            lines.append(f'                     /{key}="{value}"')

    # ORIGIN
    lines.append("ORIGIN")
    for i in range(0, n, 60):
        chunk = seq[i:i + 60]
        # 每 10 个碱基一组
        groups = [chunk[j:j + 10] for j in range(0, len(chunk), 10)]
        line_num = f"{i + 1:>9}"
        lines.append(f"{line_num} {' '.join(groups)}")

    lines.append("//")

    genbank_text = "\n".join(lines)

    return {
        "format": "GenBank",
        "name": name,
        "sequence_length": n,
        "features_count": len(features) + 1,  # +1 for source
        "content": genbank_text,
    }


def export_fragment_list(fragments: list, format: str = "text") -> dict:
    """
    导出片段清单。

    Args:
        fragments: 片段列表
        format: "text" / "csv" / "json"
    """
    if format == "csv":
        header = "Name,Start,End,Length,GC%,Overlap_5_Tm,Overlap_3_Tm"
        rows = [header]
        for f in fragments:
            rows.append(
                f"{f.get('name','')},{f.get('start','')},{f.get('end','')},{f.get('length','')},"
                f"{f.get('gc_content','')},{f.get('overlap_5_Tm_C','')},{f.get('overlap_3_Tm_C','')}"
            )
        content = "\n".join(rows)

    elif format == "json":
        import json
        content = json.dumps(fragments, indent=2)

    else:  # text
        lines = ["=" * 70, "FRAGMENT ASSEMBLY LIST", "=" * 70, ""]
        for f in fragments:
            lines.append(f"Fragment: {f.get('name', 'N/A')}")
            lines.append(f"  Position: {f.get('start', '?')}-{f.get('end', '?')}")
            lines.append(f"  Length: {f.get('length', '?')} bp")
            lines.append(f"  GC Content: {f.get('gc_content', 'N/A')}")
            ol5_tm = f.get('overlap_5_Tm_C')
            ol3_tm = f.get('overlap_3_Tm_C')
            if ol5_tm:
                lines.append(f"  5' Overlap Tm: {ol5_tm}°C")
            if ol3_tm:
                lines.append(f"  3' Overlap Tm: {ol3_tm}°C")
            ol3 = f.get('overlap_3prime', '')
            if ol3:
                lines.append(f"  3' Overlap Seq: {ol3}")
            lines.append("")
        lines.append("=" * 70)
        content = "\n".join(lines)

    return {
        "format": format,
        "fragment_count": len(fragments),
        "content": content,
    }
