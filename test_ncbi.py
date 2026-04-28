import sys
sys.path.insert(0, r"d:\gener")
from tools.ncbi_search import search_nucleotide
r = search_nucleotide("CEN3 centromere", max_results=3, organism="Saccharomyces cerevisiae")
print("Count:", r.get("total_count"))
for x in r.get("results", []):
    print(f"  {x['accession']} | {x['length']}bp | {x['title'][:80]}")
