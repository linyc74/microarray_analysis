import os
import pandas as pd
from os.path import abspath
from typing import List, Any
from .template import Processor


class GSEA(Processor):

    intensity_df: pd.DataFrame
    gene_info_df: pd.DataFrame
    sample_info_df: pd.DataFrame
    gene_name_column: str
    sample_group_column: str
    control_group_name: str
    experimental_group_name: str
    gene_sets_gmt: str

    expression_txt: str
    groups_cls: str

    def main(
            self,
            intensity_df: pd.DataFrame,
            gene_info_df: pd.DataFrame,
            sample_info_df: pd.DataFrame,
            gene_name_column: str,
            sample_group_column: str,
            control_group_name: str,
            experimental_group_name: str,
            gene_sets_gmt: str):

        self.intensity_df = intensity_df
        self.gene_info_df = gene_info_df
        self.sample_info_df = sample_info_df
        self.gene_name_column = gene_name_column
        self.sample_group_column = sample_group_column
        self.control_group_name = control_group_name
        self.experimental_group_name = experimental_group_name
        self.gene_sets_gmt = gene_sets_gmt

        self.build_expression_txt()
        self.build_groups_cls()
        self.run_gsea()

    def build_expression_txt(self):
        self.expression_txt = BuildExpressionTxt(self.settings).main(
            intensity_df=self.intensity_df,
            gene_info_df=self.gene_info_df,
            gene_name_column=self.gene_name_column)

    def build_groups_cls(self):
        self.groups_cls = BuildGroupsCls(self.settings).main(
            intensity_df=self.intensity_df,
            sample_info_df=self.sample_info_df,
            sample_group_column=self.sample_group_column,
            control_group_name=self.control_group_name,
            experimental_group_name=self.experimental_group_name)

    def run_gsea(self):
        RunGSEA(self.settings).main(
            expression_txt=self.expression_txt,
            groups_cls=self.groups_cls,
            gene_sets_gmt=self.gene_sets_gmt,
            control_group_name=self.control_group_name,
            experimental_group_name=self.experimental_group_name,
            out_dirname='gsea'
        )


class BuildExpressionTxt(Processor):

    intensity_df: pd.DataFrame
    gene_info_df: pd.DataFrame
    gene_name_column: str

    output_txt: str

    def main(
            self,
            intensity_df: pd.DataFrame,
            gene_info_df: pd.DataFrame,
            gene_name_column: str) -> str:

        self.intensity_df = intensity_df.copy()
        self.gene_info_df = gene_info_df.copy()
        self.gene_name_column = gene_name_column

        self.merge_gene_info()
        self.drop_genes_without_name()
        self.set_gene_name_as_index()
        self.add_empty_description_column()
        self.write_expression_txt()

        return self.output_txt

    def merge_gene_info(self):
        self.intensity_df = self.intensity_df.merge(
            right=self.gene_info_df[[self.gene_name_column]],
            right_index=True,
            left_index=True,
            how='left'
        )

    def drop_genes_without_name(self):
        n = len(self.intensity_df)
        self.intensity_df = self.intensity_df[self.intensity_df[self.gene_name_column].notna()]
        msg = f'For GSEA, drop genes without name (i.e. symbol), {n} -> {len(self.intensity_df)}'
        self.logger.info(msg)

    def set_gene_name_as_index(self):
        self.intensity_df = self.intensity_df.set_index(self.gene_name_column, drop=True)
        self.intensity_df.index.name = 'Name'

    def add_empty_description_column(self):
        self.intensity_df['Description'] = 'na'
        columns = self.intensity_df.columns.to_list()
        reordered = columns[-1:] + columns[0:-1]
        self.intensity_df = self.intensity_df[reordered]

    def write_expression_txt(self):
        self.output_txt = f'{self.workdir}/gsea-expression.txt'
        self.intensity_df.to_csv(self.output_txt, sep='\t', index=True)


class BuildGroupsCls(Processor):

    intensity_df: pd.DataFrame
    sample_info_df: pd.DataFrame
    sample_group_column: str
    control_group_name: str
    experimental_group_name: str

    n_total_samples: int
    n_unique_groups: int
    sample_group_names: List[str]
    cls_text: str
    output_cls: str

    def main(
            self,
            intensity_df: pd.DataFrame,
            sample_info_df: pd.DataFrame,
            sample_group_column: str,
            control_group_name: str,
            experimental_group_name: str) -> str:

        self.intensity_df = intensity_df
        self.sample_info_df = sample_info_df
        self.sample_group_column = sample_group_column
        self.control_group_name = control_group_name
        self.experimental_group_name = experimental_group_name

        self.set_numbers()
        self.set_sample_group_names()
        self.set_cls_text()
        self.write_output_cls()

        return self.output_cls

    def set_numbers(self):
        self.n_total_samples = len(self.intensity_df.columns)
        self.n_unique_groups = len(self.sample_info_df[self.sample_group_column].unique())

    def set_sample_group_names(self):
        self.sample_group_names = []
        for sample in self.intensity_df.columns:
            group = self.sample_info_df.loc[sample, self.sample_group_column]
            self.sample_group_names.append(group)

    def set_cls_text(self):
        a = ' '.join(unique(self.sample_group_names))
        b = ' '.join(self.sample_group_names)
        self.cls_text = f'''\
{self.n_total_samples} {self.n_unique_groups} 1
# {a}
{b}'''

    def write_output_cls(self):
        self.output_cls = f'{self.workdir}/gsea-groups.cls'
        with open(self.output_cls, 'w') as fh:
            fh.write(self.cls_text)


def unique(lst: List[Any]) -> List[Any]:
    ret = []
    for item in lst:
        if item not in ret:
            ret.append(item)
    return ret


class RunGSEA(Processor):

    ANALYSIS_NAME = 'gsea'
    ENRICHMENT_STATISTIC = 'weighted'
    METRIC_FOR_RANKING_GENES = 'Signal2Noise'
    GENE_LIST_SORTING_MODE = 'real'
    GENE_LIST_ORDERING_MODE = 'descending'
    COLLAPSE_REMAP_TO_GENE_SYMBOLS = 'No_Collapse'
    COLLAPSING_MODE_FOR_PROBE_SETS_GREATER_OR_EQUAL_THAN_1_GENE = 'Max_probe'
    NORMALIZATION_MODE = 'meandiv'
    NUMBER_OF_PERMUTATIONS = 1000
    PERMUTATION_TYPE = 'phenotype'
    SEED_FOR_PERMUTATION = 149
    RANDOMIZATION_MODE = 'no_balance'
    CREATE_GCT_FILES = 'false'
    CREATE_SVG_PLOT_IMAGES = 'false'
    OMIT_FEATURES_WITH_NO_SYMBOL_MATCH = 'true'
    MAKE_DETAILED_GENE_SET_REPORT = 'true'
    MEDIAN_FOR_CLASS_METRICS = 'false'
    NUMBER_OF_MARKERS = 100
    PLOT_GRAPHS_FOR_THE_TOP_SETS_OF_EACH_PHENOTYPE = 20
    SAVE_RANDOM_RANKED_LISTS = 'false'
    MAX_SIZE_EXCLUDE_LARGER_SETS = 500
    MIN_SIZE_EXCLUDE_SMALLER_SETS = 15
    MAKE_A_ZIPPED_FILE_WITH_ALL_REPORTS = 'false'

    expression_txt: str
    groups_cls: str
    gene_sets_gmt: str
    control_group_name: str
    experimental_group_name: str
    out_dirname: str

    args: List[str]

    def main(
            self,
            expression_txt: str,
            groups_cls: str,
            gene_sets_gmt: str,
            control_group_name: str,
            experimental_group_name: str,
            out_dirname: str):

        self.expression_txt = expression_txt
        self.groups_cls = groups_cls
        self.gene_sets_gmt = gene_sets_gmt
        self.control_group_name = control_group_name
        self.experimental_group_name = experimental_group_name
        self.out_dirname = out_dirname

        self.make_all_paths_absolute()
        self.set_args()
        self.run_gsea()

    def make_all_paths_absolute(self):
        self.expression_txt = abspath(self.expression_txt)
        self.groups_cls = abspath(self.groups_cls)
        self.gene_sets_gmt = abspath(self.gene_sets_gmt)
        self.workdir = abspath(self.workdir)
        self.outdir = abspath(self.outdir)

    def set_args(self):
        self.args = [
            'gsea-cli.sh GSEA',
            f'-res {self.expression_txt}',
            f'-cls {self.groups_cls}#{self.experimental_group_name}_versus_{self.control_group_name}',
            f'-gmx {self.gene_sets_gmt}',
            f'-out {self.outdir}/{self.out_dirname}',
            f'-collapse {self.COLLAPSE_REMAP_TO_GENE_SYMBOLS}',
            f'-mode {self.COLLAPSING_MODE_FOR_PROBE_SETS_GREATER_OR_EQUAL_THAN_1_GENE}',
            f'-norm {self.NORMALIZATION_MODE}',
            f'-nperm {self.NUMBER_OF_PERMUTATIONS}',
            f'-permute {self.PERMUTATION_TYPE}',
            f'-rnd_seed {self.SEED_FOR_PERMUTATION}',
            f'-rnd_type {self.RANDOMIZATION_MODE}',
            f'-scoring_scheme {self.ENRICHMENT_STATISTIC}',
            f'-rpt_label {self.ANALYSIS_NAME}',
            f'-metric {self.METRIC_FOR_RANKING_GENES}',
            f'-sort {self.GENE_LIST_SORTING_MODE}',
            f'-order {self.GENE_LIST_ORDERING_MODE}',
            f'-create_gcts {self.CREATE_GCT_FILES}',
            f'-create_svgs {self.CREATE_SVG_PLOT_IMAGES}',
            f'-include_only_symbols {self.OMIT_FEATURES_WITH_NO_SYMBOL_MATCH}',
            f'-make_sets {self.MAKE_DETAILED_GENE_SET_REPORT}',
            f'-median {self.MEDIAN_FOR_CLASS_METRICS}',
            f'-num {self.NUMBER_OF_MARKERS}',
            f'-plot_top_x {self.PLOT_GRAPHS_FOR_THE_TOP_SETS_OF_EACH_PHENOTYPE}',
            f'-save_rnd_lists {self.SAVE_RANDOM_RANKED_LISTS}',
            f'-set_max {self.MAX_SIZE_EXCLUDE_LARGER_SETS}',
            f'-set_min {self.MIN_SIZE_EXCLUDE_SMALLER_SETS}',
            f'-zip_report {self.MAKE_A_ZIPPED_FILE_WITH_ALL_REPORTS}',
            f'1> {os.path.abspath(self.outdir)}/gsea.log',
            f'2> {os.path.abspath(self.outdir)}/gsea.log',
        ]

    def run_gsea(self):
        cwd = os.getcwd()
        os.chdir(self.workdir)  # to make the gsea temp directory appear in workdir
        self.call(self.CMD_LINEBREAK.join(self.args))
        os.chdir(cwd)  # change back to the original directory
