# AutoSGI: Automated Feature Selection for Subgroup Identification

**Dharani et al.** AutoSGI: automated feature selection for subgroup identification (2026), *Submitted*. [link]()

---

## Problem Definition

High-dimensional omics data often contain heterogeneous and partially redundant signals that obscure meaningful subgroup structure and reduce the effectiveness of subgroup identification. 

## Proposed Solution

AutoSGI introduces the **first unified framework** for multi-scale feature selection and subgroup identification, enabling the discovery of clinically relevant biomarkers *and* patient subgroups.

## What is New? 

1. Partitioning features into compact subsets (either predefined by the user or via data-driven clustering)
2. Performing hierarchical clustering on samples within each subset
3. Testing subgroup splits against clinical outcomes
4. Correcting for multiple testing across all evaluated subgroup structures
5. Possibility to create meta-features (e.g., ssGSEA, DPT)

This enables robust identification of clinically meaningful subgroups at many resolutions. 

## Key Features

- Automated feature selection (data-driven or pathway-based)
- Hierarchical subgroup identification using SGI (https://github.com/krumsieklab/sgi)
- Association testing across multiple phenotype types
- Support for ordinal, categorical, and continuous outcomes
- Built-in visualization (tree plots, overview plots, feature tables)
- Exportable results (PDFs, Excel tables)
- Compatibility with multi-omics and meta-features (e.g., ssGSEA, DPT)

---

## Installation

Install from GitHub:

```r
install.packages("devtools")
devtools::install_github("richabatra/AutoSGI", subdir = "autosgi")
```

## Getting Started

### General workflow

AutoSGI follows a simple two-step workflow:

1. **Initialize analysis parameters** using `sgi_params_init()`
2. **Run a selection method** (`hierarchical_selection()` or `set_selection()`)

The `autosgi_params` object serves as the central container for the analysis. It stores:

- the input `dataset` (samples × features)
- clinical/phenotype data (`clins`)
- minimum cluster size and other parameters
- user-defined statistical tests (optional)

A separate `rule` object (created with `rule_init()`) defines how subgroup results are filtered and selected.

Once these objects are defined, the selection functions (`hierarchical_selection()` or `set_selection()`) execute the full workflow:

- feature selection (user-defined via `set_selection()` or data-driven via hierarchical clustering in `hierarchical_selection()`)
- construction of multiple SGI trees, one *per feature subset* (each tree contains subgroups at different hierarchical levels)
- association testing to identify all clinically distinct subgroups 
- multiple testing correction across all trees and subgroup splits
- result filtering and visualization

---

### Basic workflow 

An example with data-driven hierarchical selection:
```r
autosgi_params <- sgi_params_init(
  dataset = as.data.frame(dataset),
  clins = clins,
  minsize = nrow(dataset) / 20
)

rule <- rule_init() #default parameters

results <- hierarchical_selection(
  rule = rule,
  autosgi_params,
  cluster_min = 2,
  plot = TRUE,
  supp_plot = TRUE,
  summary_plot = TRUE,
  correction_opt = "simes"
)

Using predefined feature sets (e.g., pathways):
```r
results <- set_selection(
  rule = rule_init(),
  autosgi_params,
  feature_sets = feature_sets
)
```

## Applications
Two case studies are described in the AutoSGI manuscript and supplement. The code for the ROS/MAP and ADNI examples are in this repository to showcase example workflows and illustrate how AutoSGI can be used for clinically relevant subgroup identification (disease staging, trajectory analysis, etc).

