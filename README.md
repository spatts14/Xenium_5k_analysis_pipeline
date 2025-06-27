# Reproducible Spatial Transcriptomics Pipeline with RSE Best Practices

<!-- A brief description of your exemplar, which may include an image -->
This is a brief abstract of my exemplar, which includes a representative image.
![ST Lung cancer FFPE](docs/readme_image.pdf)

<!-- Author information -->
This exemplar was developed at Imperial College London by Sara Patti in
collaboration with Adrian D'Alessandro from Research Software Engineering and
Jesus Urtasun from Research Computing & Data Science at the Early Career
Researcher Institute.

<!-- Learning Outcomes.
Aim for 3 - 4 points that illustrate what knowledge and
skills will be gained by studying your ReCoDE exemplar. -->
## Learning Outcomes 🎓

<!-- TODO: NEED TO BE REVISED -->

After completing this exemplar, students will:

- Analyze spatial transcriptomic data and perform spatial statistics (Xenium)
- Develop a reporducible pipeline
- Implement RSE best practices (e.g testing, continuous integration)

<!-- Audience. Think broadly as to who will benefit. -->
## Target Audience 🎯

1) Biologists interested in developing bioinformatic pipelines
2) RSE interested in analyzing spatial transcriptomics data

<!-- Requirements.
What skills and knowledge will students need before starting?
- Python
- NextFlow (?)

Is it a prerequisite skill or learning outcome?
e.g. If your project uses a niche library, you could either set it as a
requirement or make it a learning outcome above. If a learning outcome,
you must include a relevant section that helps with learning this library.
-->

## Prerequisites ✅

### Academic 📚

- Basic understanding of Python programming
- Familiarity with [scverse ecosystem](https://scverse.org/) (e.g. scanpy, squidpy)
- Familiarity with data analysis and statistics (including spatial statistics)
- Familiarity with spatial transcriptomics

### System 💻

<!-- TODO: NEED TO BE REVISED -->

- System requirements (e.g. Python 3.11+, Anaconda, 50 GB disk space, etc.)
- Hardware or HPC requirements (if any)

<!-- Quick Start Guide. Tell learners how to engage with the exemplar. -->
## Getting Started 🚀

<!-- TODO: NEED TO BE REVISED -->

1. Start by downloading the [Xenium Lung FFPE data](https://www.10xgenomics.com/datasets/ffpe-human-lung-cancer-data-with-human-immuno-oncology-profiling-panel-and-custom-add-on-1-standard)
   - Data can be downloaded from the 10x Genomics website, or directly from the command line.
      - If downloading from the website, download the `Xenium_V1_Human_Lung_Cancer_Addon_FFPE_outs.zip` file.
      - If downloading from the command line, use the following command:

        ```bash
        curl -O https://cf.10xgenomics.com/samples/xenium/2.0.0/Xenium_V1_Human_Lung_Cancer_Addon_FFPE/Xenium_V1_Human_Lung_Cancer_Addon_FFPE_outs.zip
        ```

   - Unzip the downloaded file.
     - If you downloaded the file from the website, unzip it using your preferred method.
     - If you downloaded the file from the command line, use the following command:

        ```bash
        unzip Xenium_V1_Human_Lung_Cancer_Addon_FFPE_outs.zip
        ```

2. Create new virtual environment using `conda` or `venv`:
3. Run pipeline

Briefly describe how this project fits in your discipline, why you chose
to work on it, and what other disciplines may find it useful.

<!-- Software. What languages, libraries, software you use. -->
## Software Tools 🛠️

- Python
- squidpy
- MuSpAn

<!-- Repository structure. Explain how your code is structured. -->
## Project Structure 🗂️

<!-- TODO: NEED TO BE REVISED -->

Overview of code organisation and structure.

```text
.
├── notebooks
│ ├── ex1.ipynb
├── src
│ ├── file1.py
│ ├── file2.cpp
│ ├── ...
│ └── data
├── docs
└── test
```

Code is organised into logical components:

- `notebooks` for tutorials and exercises
- `src` for core code,
- `data` contains needed datasets
- `docs` for documentation
- `test` for testing scripts

## Roadmap 🗺️

### Preprocessing & Quality Control

<!-- TODO: NEED TO BE REVISED -->

Goal: Ensure clean, usable spatial gene expression data.

- Calculate quality metrics
- Filter low-quality genes and cells
- Normalize and transform gene counts

### Dimensionality Reduction & Clustering

Goal: Identify patterns and groups of similar gene expression profiles.

- Compute PCA and neighbors
- Compute and plot UMAP
- Cluster cells using Leiden algorithms
- Visualize clusters on UMAP

### Annotation & Cell Type Identification

Goal: Assign biological meaning to clusters.

- Compute differentially expressed genes for each cluster
- Visualize cluster marker genes
- Identify cell types with marker genes
- Annotate clusters with known cell types

### Spatial Mapping & Visualization

Goal: Map gene expression and clusters back to their spatial context.

- Overlay expression and clusters on tissue image
- Plot spatially enriched genes
- Map cell types or states in space

### Spatial Statistics & Spatially Variable Genes

Goal: Quantify spatial patterns and variability utilizing spatial statistics.

We use two different approaches to spatial statistics: Squidpy and MuSpAn.

- Compute spatial autocorrelation (e.g. Moran's I)

### TODO: Differential Expression & Functional Analysis

Goal: Discover meaningful biology.

- Spatially variable genes (SVGs)
- DE between regions or conditions
- Pathway or GO enrichment

<!-- Data availability (remove this section if no data used) -->
## Data 📊

<!-- TODO: NEED TO BE REVISED -->

- Small toy dataset for testing and development (TBD)
- [Xenium Lung FFPE data](https://www.10xgenomics.com/datasets/ffpe-human-lung-cancer-data-with-human-immuno-oncology-profiling-panel-and-custom-add-on-1-standard)

<!-- Best practice notes. -->
## Best Practice Notes 📝

<!-- TODO: NEED TO BE REVISED -->

- Git version control
- Virtual environments (e.g. conda, venv)
- Code modularity (e.g. functions, classes)
- Code documentation (e.g. docstrings, comments)
- Code style (e.g. PEP 8 for Python)
- Code testing
- Use of continuous integration (pre-commit, ruff) (?)

<!-- Estimate the time it will take for a learner to progress through the exemplar. -->
## Estimated Time ⏳

<!-- TODO: NEED TO BE REVISED -->

| Task       | Time    |
| ---------- | ------- |
| Reading    | 3 hours |
| Practising | 3 hours |

<!-- Any references, or other resources. -->
## Additional Resources 🔗

### Learn more about spatial transcriptomics

- [An introduction to spatial transcriptomics for biomedical research](https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-022-01075-1)
- [10x Genomics Xenium documentation](https://www.10xgenomics.com/xenium)
- [Single Cell Spatial Transcriptomics: 10x Genomics Xenium](https://www.youtube.com/watch?v=R4ppUVWjm7s)
- [Best practices for single cell  and spatial transcriptomics](https://www.sc-best-practices.org/preamble.html)

### Learn more about networks and spatial statistics

- [Network Science](https://networksciencebook.com/)
- [Graph Representation](https://www.youtube.com/watch?v=k1wraWzqtvQ)
- [Theory of Spatial statistics](https://www.paulamoraga.com/book-spatial/index.html)

### Learn more about our tools and libraries

- [30-days-of-python](https://github.com/Asabeneh/30-Days-Of-Python/tree/master)
- [scverse documentation](https://scverse.org/)
- [squidpy documentation](https://squidpy.readthedocs.io/en/stable/)
- [scanpy documentation](https://scanpy.readthedocs.io/en/stable/)
- [MuSpAn documentation](https://www.muspan.co.uk/resources)

### Video Tutorials

We have included bioinformatic bloggers that can help you get started with understanding key concepts in bioinformatics and transcriptomics analysis:

- [Sanbomics](https://www.youtube.com/@sanbomics)
- [Biostatsquid](https://www.youtube.com/@biostatsquid)
- [Bioinformagician](https://www.youtube.com/@Bioinformagician)

<!-- LICENCE.
Imperial prefers BSD-3. Please update the LICENSE.md file with the current year.
-->
## Licence 📄

This project is licensed under the [BSD-3-Clause license](https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics/blob/main/LICENSE.md).
