# Herb-marker

## Function Overview

This work utilizes the **Toolbox**, a publicly available repository that encapsulates a variety of functions for network-based analyses. Within this framework, several customized Python scripts and class modules were employed to analyze transcription factor–target (TF–T) and protein–protein interaction (PPI) data sets in the context of cancer and traditional Chinese medicine (TCM) associations.  

- **`TFT.py`** — Defines a Python class for processing and managing TF–gene network data, including data loading, network construction, and feature extraction.  
- **`PPI.py`** — Defines a Python class for handling PPI network data, supporting interaction parsing, network representation, and topological analysis.  
- **`cal_cancer_score_TFT_sim.py`** and **`cal_cancer_score_PPI_sim.py`** — Calculate association scores between **cancer and related TCM compounds** based on **similarity metrics**, quantifying how closely network features (TFs, targets, or proteins) resemble each other across disease and TCM contexts.  
- **`cal_cancer_score_TFT_dist.py`** and **`cal_cancer_score_PPI_dist.py`** — Compute association scores between **cancer and related TCM compounds** using **distance-based measures**, evaluating the **network distance** between TF–gene or protein pairs involved in cancer and TCM pathways.  

Together, these scripts and modules, combined with the analytical capabilities of the Toolbox, provide a comprehensive framework for investigating network-level associations between cancer mechanisms and traditional Chinese medicine.

---

## Research Context

The analytical workflow aims to explore the **molecular relationships between cancer mechanisms and TCM compounds** through a network-based perspective.  

1. **Data Preparation**  
   - Cancer-related genes and proteins are collected from curated biomedical databases.  
   - TCM compound targets are extracted from public resources such as TCMSP and ETCM.  
   - Both datasets are standardized and mapped to unified gene identifiers.  

2. **Network Construction**  
   - The `TFT.py` module constructs transcription factor–target (TF–T) networks for both cancer and TCM contexts.  
   - The `PPI.py` module builds protein–protein interaction (PPI) networks, integrating disease and compound targets into a unified network structure.  

3. **Association Analysis**  
   - The `cal_cancer_score_TFT_sim.py` and `cal_cancer_score_PPI_sim.py` scripts compute **similarity-based association scores** to quantify functional resemblance between cancer and TCM targets.  
   - The `cal_cancer_score_TFT_dist.py` and `cal_cancer_score_PPI_dist.py` scripts calculate **distance-based association scores**, reflecting the topological proximity of targets in the integrated networks.  

4. **Interpretation and Validation**  
   - High-scoring associations suggest potential **mechanistic links between TCM compounds and cancer pathways**, providing candidates for further biological validation or pharmacological study.  

This integrated framework facilitates systematic discovery of **cancer–TCM molecular correlations**, supporting evidence-based modernization of traditional Chinese medicine through computational network analysis.


## Citation

* If you use biomedical data base parsers or proximity related methods please cite: Guney E, Menche J, Vidal M, Barab&aacute;si AL. Network-based in silico drug efficacy screening. Nat. Commun. 7:10331 doi: 10.1038/ncomms10331 (2016). [link](http://www.nature.com/ncomms/2016/160201/ncomms10331/full/ncomms10331.html)
