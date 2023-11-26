import pandas as pd
from typing import Optional, List
from .template import Processor
from .tools import get_temp_path, left_join


class Limma(Processor):

    normalized_intensity_df: pd.DataFrame
    sample_df: pd.DataFrame
    sample_group_column: str
    control_group_name: str
    experimental_group_name: str
    probe_df: pd.DataFrame
    gene_name_column: str
    gene_description_column: Optional[str]

    normalized_intensity_csv: str
    sample_groups: List[str]
    r_script: str
    differential_expression_csv: str
    differential_expression_df: pd.DataFrame

    def main(
            self,
            normalized_intensity_df: pd.DataFrame,
            sample_df: pd.DataFrame,
            sample_group_column: str,
            control_group_name: str,
            experimental_group_name: str,
            probe_df: pd.DataFrame,
            gene_name_column: str,
            gene_description_column: Optional[str]):

        self.normalized_intensity_df = normalized_intensity_df
        self.sample_df = sample_df
        self.sample_group_column = sample_group_column
        self.control_group_name = control_group_name
        self.experimental_group_name = experimental_group_name
        self.probe_df = probe_df
        self.gene_name_column = gene_name_column
        self.gene_description_column = gene_description_column

        self.write_limma_input_csv()
        self.set_limma_output_csv()
        self.set_sample_groups()
        self.correct_illegal_characters()
        self.set_r_script()
        self.run_r_script()
        self.read_limma_output_csv()
        self.add_gene_name_and_description()
        self.rewrite_output_csv()

    def write_limma_input_csv(self):
        self.normalized_intensity_csv = get_temp_path(
            prefix=f'{self.workdir}/limma-normalized-intensity',
            suffix='.csv')
        self.normalized_intensity_df.to_csv(self.normalized_intensity_csv, index=True)

    def set_limma_output_csv(self):
        self.differential_expression_csv = f'{self.outdir}/differential-expression.csv'

    def set_sample_groups(self):
        self.sample_groups = []
        for sample in self.normalized_intensity_df.columns:
            group = self.sample_df.loc[sample, self.sample_group_column]
            self.sample_groups.append(group)

    def correct_illegal_characters(self):
        self.control_group_name = replace_illegal_characters(self.control_group_name)
        self.experimental_group_name = replace_illegal_characters(self.experimental_group_name)
        self.sample_groups = [
            replace_illegal_characters(name) for name in self.sample_groups
        ]

    def set_r_script(self):
        sample_conditions = ', '.join([f'"{x}"' for x in self.sample_groups])
        conditions = ', '.join([f'"{x}"' for x in unique(self.sample_groups)])

        self.r_script = f'''\
library(limma)

data <- read.csv(
  "{self.normalized_intensity_csv}",
  row.names = 1
)

# group label for each sample (i.e. column)
conditions <- c({sample_conditions})

# create design matrix
design <- model.matrix(~ 0 + factor(conditions))
# rename columns of the design matrix
colnames(design) <- c({conditions})

fit <- lmFit(data, design)

contrast.matrix <- makeContrasts(
  {self.experimental_group_name} - {self.control_group_name},
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

# adjust for multiple testing using FDR, sort by the B-statistic
topTable(fit2, adjust="fdr", sort.by="B")

results <- topTable(fit2, number=Inf, adjust="fdr", sort.by="B")
write.csv(results, file="{self.differential_expression_csv}")'''

    def run_r_script(self):
        r_file = f'{self.outdir}/limma.R'
        with open(r_file, 'w') as fh:
            fh.write(self.r_script)

        log = f'{self.outdir}/limma.log'
        cmd = self.CMD_LINEBREAK.join([
            'Rscript',
            r_file,
            f'1> {log}',
            f'2> {log}'
        ])
        self.call(cmd)

    def read_limma_output_csv(self):
        self.differential_expression_df = pd.read_csv(
            self.differential_expression_csv, index_col=0)

    def add_gene_name_and_description(self):
        cols = [self.gene_name_column]
        if self.gene_description_column is not None:
            cols.append(self.gene_description_column)

        df = left_join(
            left=self.differential_expression_df,
            right=self.probe_df[cols]
        )

        columns = df.columns.tolist()
        if self.gene_description_column is not None:
            reordered = columns[-2:] + columns[:-2]
        else:
            reordered = columns[-1:] + columns[:-1]
        df = df[reordered]

        self.differential_expression_df = df

    def rewrite_output_csv(self):
        self.differential_expression_df.to_csv(
            self.differential_expression_csv, index=True)


def replace_illegal_characters(name: str) -> str:
    return name.replace(' ', '_').replace('-', '_').replace('+', '_')


def unique(list_: list) -> list:
    ret = []
    for item in list_:
        if item not in ret:
            ret.append(item)
    return ret
