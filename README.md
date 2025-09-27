# Development of Interpretable Visualizations for Time-Series Gene Expression Data  

**Master’s Thesis Project**  

---

## Overview  

This repository contains the implementation of my Master’s thesis project.  
The aim of this work was to develop both **interactive** and **compact static visualizations** of inferred gene regulatory networks (GRNs) that support exploration of gene interaction dynamics over time, as well as to develop visualizations of enrichment analysis results and clustering of time series gene expression/gene regulation data.

Key contributions of this work include:  
- An **interactive dashboard** for exploring GRNs and enrichment analyses.  
- **Static GRN visualizations** for concise representation of dynamic GRNs
- **Clustering visualizations** for comparing temporal dynamics of gene expression and regulatory interactions.  

---

## Raw Data  

All raw data used in this project is stored in the **`Data/`** directory.  

- The directory contains raw GRN data, enrichment results, and gene expression/regulation data for clustering.  
- These files **must be preprocessed** before being used in the dashboard or static visualizations.  
- Preprocessing is done via the `data_preprocessing.ipynb` notebooks (different versions exist for human GRN, yeast GRN, enrichment, and clustering datasets).  

### Preprocessed Data Formats (used as input for visualizations)

- **GRN data** (consolidated dataframe of interactions across all modules):  
  - `Regulator` *(string)* – regulator gene name  
  - `TargetGene` *(string)* – target gene name  
  - `Weight` *(float)* – interaction strength coefficient  
  - `Module` *(int or string)* – time frame/module of the interaction  
  - `RegulatorBiotype` *(string, optional)* – biotype of regulator gene  
  - `TargetBiotype` *(string, optional)* – biotype of target gene  

- **Enrichment results data**:  
  - `Gene` *(string)* – gene name  
  - `Path_ID` *(string)* – pathway identifier  
  - `Path_Description` *(string)* – biological pathway description  
  - `Module` *(int or string)* – module of origin  
  - `Source` *(string)* – database used for pathway annotation  

- **Clustering input data**:  
  - Gene expression: `clustering/input_data/Cluster_F10/Genexpression/`  
  - Gene regulation (consolidated file): `clustering/input_data/Yeast_consolidated_gene_regulation_network.csv`  

---

## Installation  

First, install all dependencies:  

```bash
pip install -r requirements.txt
```

## Components  

### 1. Interactive Dashboard  

The interactive dashboard provides **visualizations of GRNs across different modules (time frames)** and supports enrichment analysis exploration.  
The dashboard code is located in: `streamlit_dashboard/streamlit_dashboard_with_dataupload.py`

- Implemented with **Streamlit**  
- Designed with a **modular architecture** and **decision-based workflow** that adapts dynamically to user input.  
- Supports three primary scenarios:  
  1. Uploading only a GRN file → GRN view  
  2. Uploading only an enrichment file → Enrichment view  
  3. Uploading both files → Combined view  

**Views:**  
- **GRN View**: Interactive exploration of regulatory interactions between genes.  
- **Enrichment View**: Visualization of pathway enrichment results using interactive plots.  

**Run locally:**  
- run locally with command 
```bash 
streamlit run <path_to_the_dashboard>
```
for example:
```bash
streamlit run streamlit_dashboard/streamlit_dashboard_with_dataupload.py
```

### 2. Static GRN Visualizations  

The static visualizations provide a **concise overview** of network structures and hub genes. They are compact and suitable for publications.  

- Implemented using **Graphviz** layouts:  
  - Radial (`twopi`) layout  
  - Dot layout  
- Scripts are located in the **`static_graphs/`** directory.  

---

### 3. Clustering Visualizations  

Clustering results allow exploration of gene expression and regulatory dynamics across modules.  

**Input datasets:**  
- Gene expression: `clustering/input_data/Cluster_F10/Genexpression/`  
- Gene regulation: `clustering/input_data/Yeast_consolidated_gene_regulation_network.csv`  
- scripts for visualizing clustering results across modules are provided in `clustering/cluster_visualisation.ipynb`
