# Microarray Analysis

**Affymetrix microarray analysis**

## Usage

```bash
git clone https://github.com/linyc74/microarray_analysis.git
```

The cloned repository is directly executable as a python script.

```bash
python microarray_analysis \
  --normalized-intensity-table PATH/TO/INTENSITY.CSV \
  --sample-table PATH/TO/SAMPLES.CSV \
  --probe-table PATH/TO/PROBES.CSV \
  --gene-sets-gmt PATH/TO/GSEA_GENE_SETS.GMT
```

## Environment

Linux environment dependencies:
- [`GSEA`](https://www.gsea-msigdb.org/gsea/downloads.jsp)

Python:
- `pandas`
- `seaborn`
- `sklearn`

R:
- `limma`
