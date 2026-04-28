"""
NCBI 数据库检索工具 - NCBI Search
通过 Entrez API 检索 NCBI 核酸/蛋白质数据库，获取序列和注释信息
"""

import time
from typing import Optional
from Bio import Entrez, SeqIO
from io import StringIO

import os

# NCBI 要求提供邮箱
Entrez.email = os.getenv("NCBI_EMAIL", "your.email@example.com")
Entrez.tool = "Gener-DNA-Design-Agent"

# 请求间隔，遵守 NCBI 限速 (3 requests/sec without API key)
_LAST_REQUEST_TIME = 0
_MIN_INTERVAL = 0.4  # 秒


def _rate_limit():
    """NCBI 限速控制"""
    global _LAST_REQUEST_TIME
    now = time.time()
    elapsed = now - _LAST_REQUEST_TIME
    if elapsed < _MIN_INTERVAL:
        time.sleep(_MIN_INTERVAL - elapsed)
    _LAST_REQUEST_TIME = time.time()


def search_nucleotide(query: str, max_results: int = 10, organism: Optional[str] = None) -> dict:
    """
    在 NCBI Nucleotide 数据库中搜索序列。

    Args:
        query: 搜索关键词（基因名、功能描述、accession 等）
        max_results: 最大返回条数 (1-20)
        organism: 限定物种，如 "Saccharomyces cerevisiae"
    
    Returns:
        dict: 搜索结果列表，含 accession、描述、长度等
    """
    max_results = min(max(1, max_results), 20)
    
    # 构建搜索词
    search_term = query
    if organism:
        search_term = f"{query} AND {organism}[Organism]"

    _rate_limit()
    try:
        handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=max_results, sort="relevance")
        search_results = Entrez.read(handle)
        handle.close()
    except Exception as e:
        return {"error": f"NCBI search failed: {str(e)}", "query": search_term}

    id_list = search_results.get("IdList", [])
    total_count = int(search_results.get("Count", 0))

    if not id_list:
        return {
            "query": search_term,
            "total_count": total_count,
            "results": [],
            "note": "No results found. Try different keywords or broader search terms.",
        }

    # 获取摘要信息
    _rate_limit()
    try:
        handle = Entrez.esummary(db="nucleotide", id=",".join(id_list))
        summaries = Entrez.read(handle)
        handle.close()
    except Exception as e:
        return {"error": f"NCBI summary fetch failed: {str(e)}", "ids": id_list}

    results = []
    for item in summaries:
        results.append({
            "accession": item.get("AccessionVersion", item.get("Caption", "")),
            "title": item.get("Title", ""),
            "length": int(item.get("Length", 0)),
            "organism": item.get("Organism", ""),
            "gi": str(item.get("Gi", "")),
            "create_date": item.get("CreateDate", ""),
            "update_date": item.get("UpdateDate", ""),
        })

    return {
        "query": search_term,
        "total_count": total_count,
        "returned": len(results),
        "results": results,
    }


def fetch_sequence(accession: str, seq_start: Optional[int] = None, seq_end: Optional[int] = None, format: str = "fasta") -> dict:
    """
    通过 Accession Number 从 NCBI 获取序列。

    Args:
        accession: NCBI Accession Number（如 "NC_001133.9", "NM_001301302"）
        seq_start: 起始位置 (1-based)，用于截取片段
        seq_end: 终止位置 (1-based)
        format: 返回格式 "fasta" 或 "genbank"

    Returns:
        dict: 序列信息，含序列文本、长度、注释
    """
    ret_type = "fasta" if format == "fasta" else "gb"

    kwargs = {"db": "nucleotide", "id": accession, "rettype": ret_type, "retmode": "text"}
    if seq_start is not None and seq_end is not None:
        kwargs["seq_start"] = seq_start
        kwargs["seq_stop"] = seq_end

    _rate_limit()
    try:
        handle = Entrez.efetch(**kwargs)
        raw_text = handle.read()
        handle.close()
    except Exception as e:
        return {"error": f"NCBI fetch failed for {accession}: {str(e)}"}

    # 解析序列
    try:
        record = SeqIO.read(StringIO(raw_text), "fasta" if ret_type == "fasta" else "genbank")
        seq_str = str(record.seq)
        
        result = {
            "accession": accession,
            "id": record.id,
            "name": record.name,
            "description": record.description,
            "length": len(seq_str),
            "sequence": seq_str if len(seq_str) <= 50000 else seq_str[:50000] + f"... [TRUNCATED, total {len(seq_str)}bp]",
            "format": format,
        }

        # GenBank 格式额外解析 features
        if ret_type == "gb":
            features = []
            for feat in record.features[:50]:  # 最多50个 feature
                qualifiers = {}
                for key, val in feat.qualifiers.items():
                    qualifiers[key] = val[0] if len(val) == 1 else val
                features.append({
                    "type": feat.type,
                    "location": str(feat.location),
                    "qualifiers": qualifiers,
                })
            result["features"] = features
            result["features_count"] = len(record.features)
            if hasattr(record, "annotations"):
                result["organism"] = record.annotations.get("organism", "")
                result["taxonomy"] = record.annotations.get("taxonomy", [])

        return result

    except Exception as e:
        # 如果解析失败，返回原始文本
        return {
            "accession": accession,
            "raw_text": raw_text[:10000] if len(raw_text) > 10000 else raw_text,
            "parse_error": str(e),
        }


def search_gene(query: str, organism: Optional[str] = None, max_results: int = 10) -> dict:
    """
    在 NCBI Gene 数据库中搜索基因信息。

    Args:
        query: 基因名或关键词
        organism: 限定物种
        max_results: 最大返回条数

    Returns:
        dict: 基因信息列表
    """
    max_results = min(max(1, max_results), 20)

    search_term = query
    if organism:
        search_term = f"{query} AND {organism}[Organism]"

    _rate_limit()
    try:
        handle = Entrez.esearch(db="gene", term=search_term, retmax=max_results, sort="relevance")
        search_results = Entrez.read(handle)
        handle.close()
    except Exception as e:
        return {"error": f"NCBI Gene search failed: {str(e)}"}

    id_list = search_results.get("IdList", [])
    total_count = int(search_results.get("Count", 0))

    if not id_list:
        return {"query": search_term, "total_count": 0, "results": []}

    # 获取基因摘要
    _rate_limit()
    try:
        handle = Entrez.esummary(db="gene", id=",".join(id_list))
        summaries = Entrez.read(handle)
        handle.close()
    except Exception as e:
        return {"error": f"NCBI Gene summary failed: {str(e)}"}

    results = []
    doc_sums = summaries.get("DocumentSummarySet", {}).get("DocumentSummary", [])
    for item in doc_sums:
        results.append({
            "gene_id": item.attributes.get("uid", "") if hasattr(item, "attributes") else "",
            "name": item.get("Name", ""),
            "description": item.get("Description", ""),
            "organism": item.get("Organism", {}).get("ScientificName", "") if isinstance(item.get("Organism"), dict) else str(item.get("Organism", "")),
            "chromosome": item.get("Chromosome", ""),
            "map_location": item.get("MapLocation", ""),
            "gene_type": item.get("NomenclatureSymbol", ""),
            "summary": item.get("Summary", "")[:500],
        })

    return {
        "query": search_term,
        "total_count": total_count,
        "returned": len(results),
        "results": results,
    }


def fetch_gene_sequence(gene_id: str) -> dict:
    """
    通过 Gene ID 获取基因的核酸序列（链接到 Nucleotide 数据库）。

    Args:
        gene_id: NCBI Gene ID

    Returns:
        dict: 关联的核酸序列 accession 列表
    """
    _rate_limit()
    try:
        # 通过 elink 找到关联的核酸序列
        handle = Entrez.elink(dbfrom="gene", db="nucleotide", id=gene_id, linkname="gene_nuccore_refseqrna")
        link_results = Entrez.read(handle)
        handle.close()
    except Exception as e:
        return {"error": f"NCBI elink failed: {str(e)}"}

    linked_ids = []
    for linkset in link_results:
        for link_db in linkset.get("LinkSetDb", []):
            for link in link_db.get("Link", []):
                linked_ids.append(link["Id"])

    if not linked_ids:
        # 尝试 genomic link
        _rate_limit()
        try:
            handle = Entrez.elink(dbfrom="gene", db="nucleotide", id=gene_id, linkname="gene_nuccore")
            link_results = Entrez.read(handle)
            handle.close()
            for linkset in link_results:
                for link_db in linkset.get("LinkSetDb", []):
                    for link in link_db.get("Link", [])[:5]:
                        linked_ids.append(link["Id"])
        except Exception:
            pass

    if not linked_ids:
        return {"gene_id": gene_id, "linked_sequences": [], "note": "No linked nucleotide sequences found."}

    # 获取这些序列的摘要
    _rate_limit()
    try:
        handle = Entrez.esummary(db="nucleotide", id=",".join(linked_ids[:10]))
        summaries = Entrez.read(handle)
        handle.close()
    except Exception as e:
        return {"error": f"NCBI summary failed: {str(e)}", "linked_ids": linked_ids}

    sequences = []
    for item in summaries:
        sequences.append({
            "accession": item.get("AccessionVersion", item.get("Caption", "")),
            "title": item.get("Title", ""),
            "length": int(item.get("Length", 0)),
        })

    return {
        "gene_id": gene_id,
        "linked_sequences": sequences,
        "total_linked": len(linked_ids),
    }


def blast_short(sequence: str, database: str = "nt", max_hits: int = 5) -> dict:
    """
    对短序列执行简化版 BLAST 比对（通过 NCBI BLAST API）。
    注意：BLAST 耗时较长（30-120秒），仅建议用于短片段(<1000bp)。

    Args:
        sequence: 查询 DNA 序列 (<1000bp)
        database: 数据库 "nt"(nucleotide) 或 "refseq_rna"
        max_hits: 最大命中数

    Returns:
        dict: BLAST 结果摘要
    """
    if len(sequence) > 1000:
        return {
            "error": "Sequence too long for online BLAST (max 1000bp). Use local BLAST for longer sequences.",
            "sequence_length": len(sequence),
        }

    if len(sequence) < 20:
        return {"error": "Sequence too short for meaningful BLAST results (min 20bp)."}

    from Bio.Blast import NCBIWWW, NCBIXML

    _rate_limit()
    try:
        handle = NCBIWWW.qblast(
            "blastn", database, sequence,
            hitlist_size=max_hits,
            expect=10,
            word_size=11,
        )
        blast_records = NCBIXML.read(handle)
        handle.close()
    except Exception as e:
        return {"error": f"BLAST failed: {str(e)}"}

    hits = []
    for alignment in blast_records.alignments[:max_hits]:
        best_hsp = alignment.hsps[0]
        hits.append({
            "title": alignment.title[:200],
            "accession": alignment.accession,
            "length": alignment.length,
            "e_value": best_hsp.expect,
            "score": best_hsp.score,
            "identities": best_hsp.identities,
            "align_length": best_hsp.align_length,
            "identity_pct": round(best_hsp.identities / best_hsp.align_length * 100, 1) if best_hsp.align_length else 0,
            "query_start": best_hsp.query_start,
            "query_end": best_hsp.query_end,
            "subject_start": best_hsp.sbjct_start,
            "subject_end": best_hsp.sbjct_end,
        })

    return {
        "query_length": len(sequence),
        "database": database,
        "total_hits": len(hits),
        "hits": hits,
        "note": "Results from NCBI BLAST. E-value < 1e-5 indicates significant homology.",
    }
