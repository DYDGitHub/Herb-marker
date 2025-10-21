# Herb-marker

This work utilizes the **Toolbox**, a publicly available repository that encapsulates a variety of functions for network-based analyses. Within this framework, several customized Python scripts and class modules were employed to analyze transcription factor–target (TF–T) and protein–protein interaction (PPI) data sets in cancer research.  

- **`TFT.py`** — Defines a Python class for processing and managing TF–gene network data, including data loading, network construction, and feature extraction.  
- **`PPI.py`** — Defines a Python class for handling PPI network data, supporting interaction parsing, network representation, and topological analysis.  
- **`cal_cancer_score_TFT_sim.py`** and **`cal_cancer_score_PPI_sim.py`** — Calculate cancer-related scores based on **similarity metrics**, quantifying how closely nodes (TFs, targets, or proteins) resemble each other in network topology or function.  
- **`cal_cancer_score_TFT_dist.py`** and **`cal_cancer_score_PPI_dist.py`** — Compute cancer-related scores using **distance-based measures**, assessing the **network distance** between TF–gene or protein pairs.  

Together, these scripts and modules, combined with the analytical capabilities of the Toolbox, provide a comprehensive framework for evaluating and comparing network characteristics to identify potential cancer-associated regulatory and interaction patterns.

## Citation

* If you use biomedical data base parsers or proximity related methods please cite: Guney E, Menche J, Vidal M, Barab&aacute;si AL. Network-based in silico drug efficacy screening. Nat. Commun. 7:10331 doi: 10.1038/ncomms10331 (2016). [link](http://www.nature.com/ncomms/2016/160201/ncomms10331/full/ncomms10331.html)
