#!/usr/bin/env python3
"""
Differential expression analysis with PyDESeq2 – CSV output only.
"""

import pandas as pd
import numpy as np
import os
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

# USER SETTINGS
counts_file = "full_gene_counts.csv"          
metadata_file = "metadata.csv"
output_dir = "."                              

# Change this contrast to your actual groups: ("condition", "treatment", "control")
contrast = ("condition", "male", "female")    

os.makedirs(output_dir, exist_ok=True)

# Read data
counts_df = pd.read_csv(counts_file, index_col=0)  
metadata = pd.read_csv(metadata_file, index_col=0)  

# Round fractional counts to integers (needed because featureCounts used --fraction)
counts_df = counts_df.round().astype(int)

# Transpose to samples × genes (required by PyDESeq2)
counts_df = counts_df.T

# Ensure sample order matches metadata
assert all(counts_df.index == metadata.index), "Sample names in counts and metadata don't match"

# Filter lowly expressed genes BEFORE creating the DESeq2 object
keep_genes = counts_df.sum(axis=0) >= 10 
counts_df = counts_df.loc[:, keep_genes]

# Build DESeqDataSet
design = "~ condition"
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design=design,
)

# Run DESeq2
dds.deseq2()

# Perform the statistical contrast
stat = DeseqStats(dds, contrast=contrast)
stat.summary()
res = stat.results_df
res = res.sort_values("padj")

# Save full results
res.to_csv(os.path.join(output_dir, "deseq2_results.csv"))

# Extract significant genes (padj < 0.05 and |log2FoldChange| > 1)
sig = res[(res.padj < 0.05) & (abs(res.log2FoldChange) > 1)]
sig.to_csv(os.path.join(output_dir, "deseq2_significant.csv"))

print("Analysis complete. Results saved in:", output_dir)