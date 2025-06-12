# Reproducible Spatial Transcriptomics Pipeline with RSE Best Practices

<!-- A brief description of your exemplar, which may include an image -->
This is a brief abstract of my exemplar, which includes a representative image.
![Scikit Camera Image](docs/assets/nextflow-logo.png)

<!-- Author information -->
This exemplar was developed at Imperial College London by Sara Patti in
collaboration with Adrian D'Alessandro from Research Software Engineering and
Jesus Urtasun from Research Computing & Data Science at the Early Career
Researcher Institute.

<!-- Learning Outcomes.
Aim for 3 - 4 points that illustrate what knowledge and
skills will be gained by studying your ReCoDE exemplar. -->
## Learning Outcomes 🎓

After completing this exemplar, students will:

- Analyze spatial transcriptomic data (Xenium)
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

- Required skills/knowledge (e.g. programming languages, libraries, theory, courses)

### System 💻

- System requirements (e.g. Python 3.11+, Anaconda, 50 GB disk space, etc.)
- Hardware or HPC requirements (if any)

<!-- Quick Start Guide. Tell learners how to engage with the exemplar. -->
## Getting Started 🚀

e.g. Step-by-step guide:

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

2. Visit the sections of this notebook in some particular order.
3. Attempt exercises `1a`, `1b`, etc.
4. Progress to advanced materials in the Github repository linked here.
5. Compare with solutions available in the `solutions` folder.

Briefly describe how this project fits in your discipline, why you chose
to work on it, and what other disciplines may find it useful.

<!-- Software. What languages, libraries, software you use. -->
## Software Tools 🛠️

Python, squidpy, MuSpAn

<!-- Repository structure. Explain how your code is structured. -->
## Project Structure 🗂️

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
- `src` for core code, potentially divided into further modules
- `data` within `src` for datasets
- `docs` for documentation
- `test` for testing scripts

## Roadmap 🗺️

### Preprocessing & Quality Control

Goal: Ensure clean, usable spatial gene expression data.

- Run xenium output through Space Ranger or Xenium tools
- Filter low-quality spots/cells
- Normalize gene counts

### Dimensionality Reduction & Clustering

Goal: Identify patterns and groups of similar gene expression profiles.

- PCA + UMAP/t-SNE
- Cluster by gene expression
- Identify cell types with marker genes

### Spatial Mapping & Visualization

Goal: Map gene expression and clusters back to their spatial context.

- Overlay expression and clusters on tissue image
- Plot spatially enriched genes
- Map cell types or states in space

### Differential Expression & Functional Analysis

Goal: Discover meaningful biology.

- Spatially variable genes (SVGs)
- DE between regions or conditions
- Pathway or GO enrichment

<!-- Data availability (remove this section if no data used) -->
## Data 📊

List datasets used with:

- Licensing info
- Where they are included (in the repo or external links)

<!-- Best practice notes. -->
## Best Practice Notes 📝

- Code testing and/or test examples
- Use of continuous integration (if any)
- Any other software development best practices

<!-- Estimate the time it will take for a learner to progress through the exemplar. -->
## Estimated Time ⏳

| Task       | Time    |
| ---------- | ------- |
| Reading    | 3 hours |
| Practising | 3 hours |

<!-- Any references, or other resources. -->
## Additional Resources 🔗

- Relevant sources, websites, images, AOB.

<!-- LICENCE.
Imperial prefers BSD-3. Please update the LICENSE.md file with the current year.
-->
## Licence 📄

This project is licensed under the [BSD-3-Clause license](https://github.com/ImperialCollegeLondon/ReCoDe-spatial-transcriptomics/blob/main/LICENSE.md).
