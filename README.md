# Temporal Transcriptomic Profiling of Alzheimer's Disease (NLGF Mouse Model)

## Project Overview
This project provides an advanced bioinformatic analysis of RNA-Seq data from the **NLGF (App-G-I) mouse model**, a second-generation model of Alzheimer's Disease. The analysis focuses on identifying gene expression trajectories across three critical life stages: **3 weeks (early), 8 weeks (onset), and 24 weeks (late-stage pathology)**.

## Key Features
- **Data Source**: NCBI GEO (GSE290305).
- **Advanced Statistical Modeling**: Utilizes **DESeq2** with a multi-factor interaction design (`~ age + sex + genotype + age:genotype`).
- **LRT & Wald Tests**: Employs Likelihood Ratio Test (LRT) to identify genes with significant temporal changes and Wald tests for point-to-point comparisons at late stages (T24).
- **Noise Reduction**: Implements **Variance Stabilizing Transformation (VST)** and **Log Fold Change Shrinkage (apeglm)** for robust visualization and biological interpretation.

## Analysis Pipeline
1. **Data Pre-processing**: Metadata harmonization and string manipulation for sample alignment.
2. **Exploratory Data Analysis**: PCA plots using VST-normalized counts to visualize genotype and age clusters.
3. **Differential Expression**: Identification of "Interaction Genes"—those whose response to the disease genotype changes significantly over time.
4. **Visualizations**: 
   - PCA showing developmental and pathological trajectories.
   - Heatmaps highlighting the top genes driving the interaction effect.

## How to Run
1. Ensure you have `R` and the `BiocManager` packages installed.
2. Place the count table and series matrix in the working directory.
3. Run the `RNA_Seq Analysis.R` script.

## Results
The analysis identifies a clear divergence in transcriptomic signatures between WT and NLGF mice that intensifies with age, particularly at the 24-week mark, highlighting key pathways involved in neuroinflammation and microglial activation.
