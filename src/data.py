from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "PanPred" / "test_data"
ANTIBIOTICS = ["CTZ", "CTX", "AMP", "AMX", "AMC", "TZP", "CXM", "CET", "GEN", "TBM", "TMP", "CIP"]


def load_data(data_dir=DATA_DIR):
    """
    Load the PanPred tables and align them on strains that appear in all three.
    Removes genes that are absent in all strains or present in all strains.
    Returns X (genes), labels (R/S), year, clusters (all SNP cutoffs).
    """
    metadata = pd.read_csv(data_dir / "Metadata.csv", index_col="Isolate")
    metadata["Year"] = metadata["Year"].astype("Int64")

    genes = pd.read_csv(data_dir / "AccessoryGene.csv", index_col=0).astype("uint8")
    population = pd.read_csv(data_dir / "PopulationStructure.csv", index_col=0)

    common_ids = sorted(set(metadata.index) & set(genes.index) & set(population.index))
    metadata, genes, population = (df.loc[common_ids] for df in (metadata, genes, population))

    assert metadata.index.equals(genes.index)
    assert metadata.index.equals(population.index)

    gene_counts = genes.sum()
    is_constant = (gene_counts == 0) | (gene_counts == len(genes))
    genes = genes.loc[:, ~is_constant]

    return genes, metadata[ANTIBIOTICS], metadata["Year"], population