"""Module size tables: one row per module, ranked, and a per-fit summary that counts giants."""
import pandas as pd

from isograph_benchmark.real_data import module_sizes as ms


def _write(root, cohort, region, method, sizes):
    d = root / cohort / region / "_m" / method
    d.mkdir(parents=True)
    genes, mods = [], []
    for k, n in enumerate(sizes):
        genes += [f"{method}_{region}_g{k}_{i}" for i in range(n)]
        mods += [f"M{k:03d}"] * n
    pd.DataFrame({"gene_id": genes, "module_id": mods}).to_parquet(d / "modules.parquet")


def test_sizes_are_ranked_and_giants_counted(tmp_path):
    _write(tmp_path, "gtex", "amygdala", "isograph_vae", [30, 950, 120])
    _write(tmp_path, "gtex", "amygdala", "wgcna_gene", [2000, 40])
    _write(tmp_path, "brainseq", "caudate", "isograph_vae", [25, 25])
    per_module, summary = ms.collect(tmp_path)

    iso = per_module[(per_module.region == "amygdala") & (per_module.method == "isograph_vae")]
    assert iso["n_genes"].tolist() == [950, 120, 30] and iso["size_rank"].tolist() == [1, 2, 3]
    assert abs(iso["frac_of_assigned"].sum() - 1.0) < 1e-3

    s = summary.set_index(["cohort", "region", "method"])
    assert s.loc[("gtex", "amygdala", "isograph_vae"), "n_modules_ge_900"] == 1
    assert s.loc[("gtex", "amygdala", "wgcna_gene"), "largest"] == 2000
    assert s.loc[("brainseq", "caudate", "isograph_vae"), "n_modules_ge_900"] == 0
    assert set(summary["method"]) == {"isograph_vae", "wgcna_gene"}  # absent methods are skipped


def test_partition_hash_ignores_row_order_but_not_assignment():
    a = pd.DataFrame({"gene_id": ["g1", "g2"], "module_id": ["M000", "M001"]})
    assert ms._partition_sha256(a) == ms._partition_sha256(a.iloc[::-1])
    assert ms._partition_sha256(a) != ms._partition_sha256(a.assign(module_id=["M001", "M000"]))
