# Benchmarking niche identification via domain segmentation for spatial transcriptomics data

## Overview

Spatial niche identification aims to partition tissue into multicellular microenvironments that are defined by coordinated cell type composition and spatially structured cellular states, using spatially resolved expression measurements. Across tissues, such microenvironments can reflect recurrent cellular neighborhoods, context-dependent state programs, and local cell–cell interactions that are not fully captured by gene expression alone or by coarse anatomical landmarks. Niches may align with anatomy in some settings, but they can also be sharply separated yet internally heterogeneous, appear as non-contiguous islands embedded within larger compartments, or vary continuously along gradients. These properties make it difficult to infer method performance from a single reference setting, motivating a benchmark that jointly probes multiple niche geometries and practical data regimes under a consistent task definition.

This repository contains the benchmarking framework and code for our paper. We evaluate 16 representative algorithms spanning probabilistic models, graph neural networks, deep generative models, and foundation models. The benchmark quantifies performance across complementary axes including agreement with reference niches (accuracy), spatial structure and boundary fidelity (connectivity), biological consistency of inferred niches (composition similarity), quality of learned embeddings (silhouette score), and computational efficiency (runtime and memory).

## Benchmark Methods

We benchmarked 16 representative algorithms categorized into four methodological families.

### Probabilistic & Statistical
*   **BayesSpace**: [Code](https://github.com/edward130603/BayesSpace) | [Paper: Spatial transcriptomics at subspot resolution with BayesSpace (Nature Biotechnology)](https://www.nature.com/articles/s41587-021-00935-2)
*   **BANKSY**: [Code](https://github.com/prabhakarlab/Banksy_py) | [Paper: BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis (Nature Genetics)](https://www.nature.com/articles/s41588-024-01664-3)
*   **MENDER**: [Code](https://github.com/yuanzhiyuan/MENDER) | [Paper: MENDER: fast and scalable tissue structure identification in spatial omics data (Nature Communications)](https://www.nature.com/articles/s41467-023-44367-9)

### GNN & Contrastive
*   **SpaGCN**: [Code](https://github.com/jianhuupenn/SpaGCN) | [Paper: SpaGCN: Integrating gene expression, spatial location and histology to identify spatial domains and spatially variable genes by graph convolutional network (Nature Methods)](https://www.nature.com/articles/s41592-021-01255-8)
*   **GraphST**: [Code](https://github.com/JinmiaoChenLab/GraphST) | [Paper: Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST (Nature Communications)](https://www.nature.com/articles/s41467-023-36796-3)
*   **STAGATE**: [Code](https://github.com/zhanglabtools/STAGATE) | [Paper: Deciphering spatial domains from spatially resolved transcriptomics with an adaptive graph attention auto-encoder (Nature Communications)](https://www.nature.com/articles/s41467-022-29439-6)
*   **CytoCommunity**: [Code](https://github.com/tanlabcode/CytoCommunity) | [Paper: Unsupervised and supervised discovery of tissue cellular neighborhoods from cell phenotypes (Nature Methods)](https://www.nature.com/articles/s41592-023-02124-2)
*   **SpaceFlow**: [Code](https://github.com/hongleir/SpaceFlow) | [Paper: Identifying multicellular spatiotemporal organization of cells with SpaceFlow (Nature Communications)](https://www.nature.com/articles/s41467-022-31739-w)

### Deep Generative
*   **NicheCompass**: [Code](https://github.com/Lotfollahi-lab/nichecompass) | [Paper: Quantitative characterization of cell niches in spatially resolved omics data (Nature Genetics)](https://www.nature.com/articles/s41588-025-02120-6)
*   **scNiche**: [Code](https://github.com/ZJUFanLab/scNiche) | [Paper: Identification and characterization of cell niches in tissue from spatial omics data at single-cell resolution (Nature Communications)](https://www.nature.com/articles/s41467-025-57029-9)
*   **SEDR**: [Code](https://github.com/JinmiaoChenLab/SEDR) | [Paper: Unsupervised spatially embedded deep representation of spatial transcriptomics (Genome Medicine)](https://link.springer.com/article/10.1186/s13073-024-01283-x)
*   **STACI**: [Code](https://github.com/uhlerlab/STACI) | [Paper: Graph-based autoencoder integrates spatial transcriptomics with chromatin images and identifies joint biomarkers for Alzheimer’s disease (Nature Communications)](https://www.nature.com/articles/s41467-022-35233-1)
*   **DeepLinc**: [Code](https://github.com/xryanglab/DeepLinc) | [Paper: De novo reconstruction of cell interaction landscapes from single-cell spatial transcriptome data with DeepLinc (Genome Biology)](https://link.springer.com/article/10.1186/s13059-022-02692-0)
*   **CellCharter**: [Code](https://github.com/CSOgroup/CellCharter) | [Paper: CellCharter reveals spatial cell niches associated with tissue remodeling and cell plasticity (Nature Genetics)](https://www.nature.com/articles/s41588-023-01588-4)

### Foundation Models
*   **Novae**: [Code](https://github.com/MICS-Lab/Novae) | [Paper: Novae: a graph-based foundation model for spatial transcriptomics data (Nature Methods)](https://www.nature.com/articles/s41592-025-02899-6)
*   **Nicheformer**: [Code](https://github.com/theislab/nicheformer) | [Paper: Nicheformer: a foundation model for single-cell and spatial omics (Nature Methods)](https://www.nature.com/articles/s41592-025-02814-z)

## Evaluation Metrics

Results are quantified using the following categorical metrics:

*   **Ground Truth Accuracy**: Measures agreement with reference annotations.
    *   ARI (Adjusted Rand Index)
    *   AMI (Adjusted Mutual Information)
    *   Homogeneity
    *   Completeness
    *   Macro-F1 Score
*   **Biological Consistency**: Assesses the biological relevance of inferred niches.
    *   Cell Type Cosine Similarity
*   **Spatial Structure**: Evaluates the spatial coherence of the partition.
    *   Spatial Connectivity
*   **Embedding Quality**: Measures the separation and compactness of the latent representation.
    *   Silhouette Score
*   **Computational Efficiency**: Assesses practical scalability.
    *   Runtime
    *   Peak Memory

## Code Structure

*   `annotation/`: Contains scripts and data for constructing the manually annotated high-resolution human lymph node reference and defining ground truths.
*   `benchmark/`: Contains the implementation and running scripts (Jupyter Notebooks) for the 16 benchmarked methods.
*   `simulation/`: Contains code for generating synthetic spatial transcriptomics data using SRTsim with controlled niche parameters.

## License

This project is released under the MIT License. See [LICENSE](LICENSE) for details.
