# Benchmarking of methods for spatial clustering on CosMx spatial transcriptomics data

**MSc Thesis in Bioinformatics (30 ECTS), Spring Semester 2026**

**Conducted under the supervision of Prof. Dyrskjøt Andersen and Tessa Jane Divita, Department of Molecular Medicine (MOMA)**

## About the repository

In this repository, you can find the computational environments used during this work, as well as the analytical pipelines I developed to analyse CosMx Spatial Transcriptomics data using the BASS, BANKSY, and SEDR methods. These pipelines were implemented in R and Python, based on the official implementations of the respective methods and adapted where necessary.

## Copyright and usage restrictions

**Copyright © 2026 Marina Ramírez Baños and Bladder Cancer group (MOMA). All rights reserved.**

The code in this repository was developed in the context of the MSc thesis *“Benchmarking of methods for spatial clustering on CosMx spatial transcriptomics data”*.

Unless otherwise indicated, the code is protected by copyright and is **not released under an open-source license**. No permission is granted to copy, modify, redistribute, republish, incorporate, or otherwise reuse the code, in whole or in part, without prior written authorization from the relevant rights holders.

**Any intended use of the code, including academic research, educational purposes, or commercial applications, requires prior written permission.**

If you wish to use or adapt any part of the code, please contact the author and the research group to describe your intended use and request authorization.

- **Author:** marina.ramirezbanos@gmail.com
- **Research group / Department:** lars@clin.au.dk

Third-party software, packages, methods, and code incorporated into or referenced by this repository remain subject to their respective copyright, licenses, and terms of use. This notice does not override any applicable third-party licenses.

## Acknowledgements

This work was conducted under the supervision of Prof. Dyrskjøt Andersen and Tessa Jane Divita at the Department of Molecular Medicine (MOMA).

## Abstract

### Background
The identification of cell-type populations and the characterization of their spatial organization within tissue are fundamental to unraveling the role of the tumor microenvironment (TME) in cancer progression. In this context, spatial transcriptomics (ST) and spatial clustering methods have emerged as powerful tools for generating tissue maps based on gene expression and cellular location. Multiple methods are available for spatially aware cell-type clustering. Nevertheless, there is a lack of benchmark studies evaluating the suitability of these methods across different data types and considering the impact of tissue architecture on their performance. Bladder cancer (BC) is a highly heterogeneous and prevalent disease, and the role of the TME remains largely understudied. Evaluating and identifying appropriate spatial clustering methods for ST data in BC could help establish methodological guidelines for future studies of TME composition and spatial organization.

### Objectives
To evaluate the performance of BANKSY, BASS, and SEDR spatial clustering methods on CosMx ST data from muscle-invasive bladder cancer (MIBC). To assess the effects of method-specific parameter settings and data preprocessing steps on performance. To identify parameter configurations that improve the performance of each method and compare them. Finally, to interpret the clustering results in the context of the intrinsic spatial architecture of the tissue.

### Methods
Tissue samples from 112 patients with MIBC were processed as formalin-fixed, paraffin-embedded tissue microarrays, yielding 298 tissue cores. Gene expression and spatial transcript locations were acquired using the CosMx Spatial Molecular Imager (Human 6K Discovery Panel). Cell segmentation was performed using membrane and nuclear protein markers. Data were preprocessed following Scanpy guidelines, and BANKSY, BASS, and SEDR were implemented according to their official tutorials. Method performance was evaluated using: (i) alignment with reference cell-type annotations (ARI, NMI, FMI, Homogeneity, Completeness, and Purity); (ii) coherence with intrinsic spatial organization (PAS, CHAOS, ASW, and SCS); (iii) robustness to data variation (sparsity, cellularity, and number of clusters); and (iv) computational efficiency (runtime and memory usage).

### Results
BANKSY exhibited better alignment with reference annotations and tissue spatial architecture than SEDR and BASS. BANKSY achieved a balance between capturing internal cluster composition and preserving structural integrity. MIBC exhibited a granular spatial structure with low global spatial continuity, and the detection of finer clusters improved alignment with reference annotations. Log-transformation of the data influenced the performance of all methods. BANKSY's performance was affected by the selection of the λ parameter, with the best performance observed at λ = 0.2. All methods were robust to variations in cellularity, whereas only BANKSY was affected by sparsity in the gene expression matrix. BANKSY showed the best runtime scalability but also the highest memory consumption.

### Conclusions
Among the methods evaluated, BANKSY was the most suitable spatial clustering method for MIBC CosMx ST data. It effectively captured the intrinsic spatial architecture of the tissue and generated fine-grained spatial clusters, resulting in improved alignment with reference annotations. The predicted clusters overlapped with niches identified in the reference annotations, suggesting that the method can detect cell-state representations enriched with spatial information. The selection of λ, gene expression matrix sparsity, and data log-transformation affected method performance. Nevertheless, BANKSY was robust to variations in cellularity and demonstrated stable runtime performance.

### Keywords
Tumor microenvironment, spatial transcriptomics, spatial clustering, bladder cancer, cell type, benchmarking.
