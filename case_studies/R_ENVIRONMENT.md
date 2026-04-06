# R Environment Documentation

## R Version
- **Minimum Required:** R >= 4.2.0
- **Tested On:** R 4.5.2 (2025-10-31)

## System Information
- **Platform:** x86_64-apple-darwin20
- **OS:** macOS Sequoia 15.7.3

## Required R Packages

### Core Analysis Packages
```r
autosgi ("0.0.1.0"")    # AutoSGI hierarchical selection methods
sgi ("0.0.0.9000")         # Supervised group identification
maplet ("1.1.2")     # Metabolomics analysis tools
```

### Data Manipulation
```r
tidyverse ("2.0.0")   # Data wrangling (includes dplyr, tidyr, ggplot2)
magrittr ("2.0.4" )   # Pipe operators (%>%, %<>%)
dplyr ("1.1.4")      # Data frame manipulation
tidyr ("1.3.2")      # Data reshaping
```

### Omics Data Structures
```r
SummarizedExperiment ("1.40.0")  # Bioconductor data container
```

### Enrichment and Pathway Analysis
```r
GSVA ("2.4.4")        # Gene Set Variation Analysis (ssGSEA)
```

### Diffusion Analysis
```r
destiny ("3.24.0")     # Diffusion maps
diffusionMap ("1.2.0") # Diffusion pseudotime
```

### Statistical Modeling
```r
rms ("8.1-0")         # Regression modeling strategies (ordinal regression)
```

### Visualization
```r
ggplot2 ("4.0.1")     # Grammar of graphics (included in tidyverse)
RColorBrewer ("1.1-3") # Color palettes
```

### File I/O
```r
openxlsx ("4.2.8.1")    # Reading/writing Excel files
readxl ("1.4.5")      # Reading Excel files
```

### String Processing
```r
stringr ("1.6.0")     # String manipulation
```

### Data Reshaping
```r
reshape ("0.8.10")    # Data reshaping (legacy)
reshape2 ("1.4.5")    # Data reshaping (melt/cast)
```

## Installation Instructions

### Step 1: Install Base R
Download and install R from [CRAN](https://cran.r-project.org/)

### Step 2: Install CRAN Packages
```r
install.packages(c(
  "tidyverse",
  "magrittr", 
  "openxlsx",
  "readxl",
  "stringr",
  "reshape",
  "reshape2",
  "RColorBrewer",
  "rms"
))
```

### Step 3: Install Bioconductor Packages
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "SummarizedExperiment",
  "GSVA",
  "destiny"
))
```

### Step 4: Install Custom Packages
```r
# Install autosgi, sgi, maplet
# Adjust based on where these packages are hosted
# Example for GitHub:
# install.packages("devtools")
# devtools::install_github("username/autosgi")
# devtools::install_github("username/sgi")
# devtools::install_github("username/maplet")
```

### Step 5: Install diffusionMap
```r
# If not available on CRAN/Bioconductor, may need manual installation
# or from archived source
```

## Script Dependencies

### 0_custom_functions.R
**Purpose:** Custom helper functions used across all analyses

**Dependencies:**
- destiny (diffusion maps)
- diffusionMap (DPT calculation)
- GSVA (pathway enrichment)
- SummarizedExperiment (data structures)
- rms (ordinal regression)
- ggplot2 (visualization)
- dplyr, tidyr, magrittr (data manipulation)
- reshape, reshape2 (data reshaping)
- sgi (SGI functions)
- openxlsx (Excel I/O)

---

### 1_rosmap-autosgi_case-study.R
**Purpose:** ROSMAP metabolomics AutoSGI analysis

**Dependencies:**
- maplet (data loading)
- sgi (SGI analysis)
- autosgi (hierarchical selection)
- tidyverse (data manipulation)
- magrittr (pipes)
- rms (statistical tests)
- 0_custom_functions.R

**Inputs:**
- rosmap-data/tmp_rosmap_brain_metabolomics_processed_medcor_data_pmicor.xlsx

**Outputs:**
- {date}_rosmap-results/rosmap-metabo-autosgi-data.rds
- {date}_rosmap-results/rosmap-metabo-all-sgi.pdf
- {date}_rosmap-results/rosmap-metabo-hierarchical-selection-*.pdf
- {date}_rosmap-results/rosmap-metabo-hierarchical-selection-*.xlsx

---

### 2_rosmap-plots.R
**Purpose:** Visualization of ROSMAP AutoSGI results

**Dependencies:**
- maplet, sgi, autosgi
- tidyverse, magrittr
- openxlsx, RColorBrewer
- rms, stringr
- 0_custom_functions.R

**Inputs:**
- Results from script 1 (rosmap-autosgi_case-study)

**Outputs:**
- Cluster-specific visualization PDFs

---

### 3_adni-autosgi_case-studies.R
**Purpose:** ADNI lipidomics AutoSGI analysis (3 approaches)

**Dependencies:**
- All packages from script 1, plus:
- GSVA (ssGSEA)
- destiny, diffusionMap (DPT)
- SummarizedExperiment (data structures)
- stringr
- 0_custom_functions.R

**Inputs:**
- adni-data/ADMCLIPIDOMICSMEIKLELABLONG_08_13_21_20Jun2024.csv
- adni-data/2024-09-26-long_metadata.xlsx
- adni-data/2024-09-26-base_metadata.xlsx

**Outputs:**
- {date}_adni-results/adni-lipidomics-autosgi-data.rds
- {date}_adni-results/adni-lipidomics-autosgi-ssgsea-data.rds
- {date}_adni-results/adni-lipidomics-autosgi-dpt-data.rds
- {date}_adni-results/adni-lipidomics-ssgsea-hierarchical-selection-*.pdf/xlsx
- {date}_adni-results/adni-lipidomics-dpt-hierarchical-selection-*.pdf/xlsx
- {date}_adni-results/adni-lipidomics-pathway-selection-*.pdf/xlsx

---

### 4_adni-plots.R
**Purpose:** Visualization of ADNI AutoSGI results

**Dependencies:**
- Same as script 2
- 0_custom_functions.R

**Inputs:**
- Results from script 3 (adni-autosgi_case-studies)

**Outputs:**
- Cluster-specific visualization PDFs for each method

---

## Script Execution Order

1. **0_custom_functions.R** - Source this in all other scripts
2. **1_rosmap-autosgi_case-study.R** - ROSMAP analysis
3. **2_rosmap-plots.R** - ROSMAP visualization
4. **3_adni-autosgi_case-studies.R** - ADNI analysis
5. **4_adni-plots.R** - ADNI visualization

Scripts 1-2 (ROSMAP) and scripts 3-4 (ADNI) are independent and can be run in parallel.

## Notes

- Custom packages (autosgi, sgi, maplet) may require specific installation instructions from package maintainers
- All scripts use `source("0_custom_functions.R")` and expect it in the same directory
- Scripts create dated output directories automatically

## Contact

For questions about package versions or installation, contact the script author or package maintainers.
