# Strip Plot Database(Python)

This Python script generates several hybrid plots by overlaying box plots and strip plots in order to group and compare different datasets. It represents the total number of significantly differentially expressed miRNAs, their relative percentage with respect to the number of miRNAs included in the analysis, and the percentage relative to the size of the database. Additionally, a complementary plot is generated to examine the effect of the methodology on the magnitude of the fold change.

---

## Input files

- **Directory containing sRNAde results (differential expression matrices)**
  - Selected by the user 
  - The folder for the assignment method must be named: 
	-	de_rcsa 
	-	de_rcadj
  - The matrices must be named: 
	-	DESeq2_[Condition 1]_vs_[Condition 2].tsv
	-	edgeR_[Condition 1]_vs_[Condition 2].tsv
	-	limma_[Condition 1]_vs_[Condition 2].tsv

- **Directory containing local copies of databases**
  - Selected by the user 
  - Each file must be named after its corresponding database.
	-	miRBase
	-	mircarta
	-	mirgenedb
---

## Requirements

- **Python 3**
- Required libraries:
  - `pandas`
  - `matplotlib`
  - `seaborn`

No additional external dependencies are required.

---

## Script usage

Run the script from the terminal or a Python environment:

```bash
python strip_plot_DB.py -i <input_file> -o <output_file> [-db <databases>] [-s <specie>]
````

Arguments
| Argument | Alternative_Argument | Required | Description | Default |
| --- | --- | --- | --- | --- |
| -i | --input | Yes | Path to the input file | - |
| -o | --output | Yes | Path to the output file | - |
| -db | --databases | No | Path to the local copies of databases | /shared/bak/TFG/serrano/db |
| -s | --specie | No | Specie | Hsa |

---

## Output

The script generates five output files:

1. **Tab-separated values file (tsv)**
	- Contains the processed data with all the values used in the graphical representations.
	
2. **Hybrid plot for significantly differentially expressed miRNAs (png)**

3. **Hybrid plot for percent of significantly differentially expressed miRNAs relative to the number of miRNAs included in the analysis (png)**

4. **Hybrid plot for percent of significantly differentially expressed miRNAs relative to the size of the database (png)**

5. **Hybrid plot to examine the effect of the methodology on the magnitude of the fold change (png)**

---

## How the script works

The script is divided into four main stages:

### 1️. File selection and loading

- Searches the user-selected folder for all differential expression matrices and stores them in a list.  

### 2. Generation of the `DataFrame`

- Each file is iterated over individually.
- The number of significant miRNAs and significantly differentially expressed miRNAs, along with their corresponding percentages, are calculated.
- The 80th percentile of the absolute log2FoldChange is calculated.
- Those values are added to a `DataFrame` containing information for all methodological combinations within a study.
- A copy of the `DataFrame` is saved as a `tsv` file.

### 3. Count the number of miRNAs in every database

- The number of miRNAs in each database is determined using the local copies of each repository.

### 4. Graphical representations

- The dataset's hybrid plots are generated.
