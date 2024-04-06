
# Placenta Project

## Introduction
This repository contains the codebase for a comprehensive study on the placenta, focusing on aspects such as cell-cell communication, differential gene expression, and pseudotime analyses among others. This project utilizes various bioinformatics tools and methods to analyze single-nucleus RNA sequencing (snRNA-seq) data for in-depth understanding of placental biology.

## Project Structure
The project is organized into several directories, each dedicated to a specific analysis or process in the study. Below is a brief overview of each directory and its contents:

### Batch Correction and UMAP Nonlinear Dimensional Reduction
- `1_prepare_sobj_x9_allcelltype_before_integrate_mnn.R`: Prepares single-cell object (sobj) for all cell types before integration using mutual nearest neighbors (MNN).
- `2_integrate_x9_allcelltype_mnn.R`: Integrates all cell types using MNN.
- `3_UMAP_x9_allcelltype.R`: Performs UMAP nonlinear dimensional reduction on the integrated cell types.

### Cell-Cell Communication Analyses
- `1_prepare_sobj_allcelltype_2GCsubtypes.R`: Prepares sobj for all cell types and two GC subtypes.
- `2_CellChat.R`: Analyzes cell-cell communication using CellChat.
- `3_generate_L-R_table.R`: Generates ligand-receptor (L-R) table for further analysis.

### Celltype Integration
- `1_prepare_sobj_x9_Wang_eLife_before_integrate_mnn.R`: Prepares sobj for integration according to Wang et al., eLife.
- `2_integrate_x9_Wang_eLife_mnn.R`: Integrates cell types based on Wang et al., eLife using MNN.
- `3_UMAP.R`: Performs UMAP on the integrated dataset.

### Data Processing of snRNA-seq Data
- `cellranger_aggr_lite.sh`: Lightweight version of Cell Ranger's `aggr` command for snRNA-seq data preprocessing.

### Identification of Differentially Expressed Genes
- `DEG_Macrophage.R`: Identifies differentially expressed genes in macrophages.

### Ontology Annotation
- `GO_GC_TB.R`: Performs Gene Ontology (GO) annotation for trophoblast cells (TB) and germ cells (GC).

### Pseudotime Analyses and Characterization of GC Subtypes
- `1_prepare_sobj_GC.R`: Prepares sobj for germ cells.
- `2_integrate_GC_mnn.R`: Integrates germ cell data using MNN.
- `3_plot_UMAP.R`: Plots UMAP for germ cell data visualization.
- `4_monocle.R`: Applies Monocle for pseudotime analysis.

### SCENIC Analysis
- Scripts for running SCENIC analysis, including data preparation, regulon generation, and heatmap updating.

## Getting Started
To use this code, clone the repository and navigate into each directory to run the scripts relevant to your analysis. Ensure you have all the necessary R packages and tools installed.

```bash
git clone <repository-url>
cd Placenta_Project
# navigate to specific directories as needed
```

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

For detailed information on each analysis, refer to the scripts within each directory. Ensure you read any associated comments or documentation provided within the scripts for a better understanding of their functionality.
