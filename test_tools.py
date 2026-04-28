"""快速功能测试"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from tools.thermodynamics import calc_tm, calc_nn_thermodynamics
from tools.sequence_validator import basic_stats, calc_gc_content
from tools.structure_analyzer import find_g_quadruplex
from tools.biosafety_checker import screen_sequence
from tools.loxp_designer import check_orthogonality
from tools.fragment_assembler import split_gibson
from tools.format_converter import to_fasta
from tools.telomere_centromere import design_telomere, design_centromere
from tools.complexity_scorer import score_complexity

def test_all():
    print("=== Test 1: Tm ===")
    r = calc_tm("ATCGATCGATCGATCGATCG")
    print(f"  Tm: {r['Tm_C']}°C  Method: {r['method']}")
    assert r["Tm_C"] > 0

    print("\n=== Test 2: NN Thermodynamics ===")
    r = calc_nn_thermodynamics("GCGATCGC")
    print(f"  dH={r['dH_cal_mol']} cal/mol  dS={r['dS_cal_mol_K']} cal/mol·K  dG37={r['dG_37_kcal_mol']} kcal/mol")
    assert r["dG_37_kcal_mol"] < 0

    print("\n=== Test 3: Basic Stats ===")
    r = basic_stats("ATCGATCG" * 100)
    print(f"  Length: {r['length']}bp  GC: {r['gc_content']}")
    assert r["length"] == 800

    print("\n=== Test 4: G-quadruplex ===")
    r = find_g_quadruplex("AAAGGGAAAGGGAAAGGGAAAGGGAAA")
    print(f"  G4 count: {r['g4_count']}  Risk: {r['overall_risk']}")
    assert r["g4_count"] >= 1

    print("\n=== Test 5: Biosafety (safe) ===")
    r = screen_sequence("ATCGATCG" * 50)
    print(f"  Safe: {r['safe']}  Alerts: {r['alerts_count']}")
    assert r["safe"]

    print("\n=== Test 6: loxP Orthogonality ===")
    r = check_orthogonality(["loxP", "lox511", "lox2272"])
    print(f"  All orthogonal: {r['all_orthogonal']}")
    assert r["all_orthogonal"]

    print("\n=== Test 7: Gibson Split ===")
    r = split_gibson("ATCGATCG" * 1000, target_size=2000)
    print(f"  Fragments: {r['total_fragments']}")
    assert r["total_fragments"] >= 3

    print("\n=== Test 8: FASTA ===")
    r = to_fasta("ATCGATCG" * 10, "test_seq")
    print(f"  Format: {r['format']}  Length: {r['sequence_length']}")
    assert ">test_seq" in r["content"]

    print("\n=== Test 9: Telomere ===")
    r = design_telomere("yeast", 500)
    print(f"  Length: {r['telomere_length']}bp  G4 risk: {r['g4_risk']}")
    assert r["telomere_length"] == 500

    print("\n=== Test 10: Centromere ===")
    r = design_centromere("yeast", "synIII")
    print(f"  Length: {r['centromere_length']}bp  AT: {r['overall_at_content']}")
    assert r["centromere_length"] > 100

    print("\n=== Test 11: Complexity Score ===")
    r = score_complexity("ATCGATCG" * 100)
    print(f"  Score: {r['complexity_score']}  Grade: {r['grade']}")
    assert 0 <= r["complexity_score"] <= 100

    print("\n" + "=" * 50)
    print("✅ All 11 tests passed!")

if __name__ == "__main__":
    test_all()
