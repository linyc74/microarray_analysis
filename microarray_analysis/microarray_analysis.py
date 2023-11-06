import os
import pandas as pd
from typing import Optional
from .pca import PCA
from .gsea import GSEA
from .limma import Limma
from .tools import get_files
from .heatmap import Heatmap
from .template import Processor


class MicroarrayAnalysis(Processor):

    normalized_intensity_table: str
    sample_table: str
    probe_table: str
    gene_sets_gmt: Optional[str]
    gene_name_column: str
    gene_description_column: Optional[str]
    heatmap_intensity_fraction: float
    sample_group_column: str
    control_group_name: str
    experimental_group_name: str

    normalized_intensity_df: pd.DataFrame
    sample_df: pd.DataFrame
    probe_df: pd.DataFrame

    def main(
            self,
            normalized_intensity_table: str,
            sample_table: str,
            probe_table: str,
            gene_sets_gmt: Optional[str],
            gene_name_column: str,
            gene_description_column: Optional[str],
            heatmap_intensity_fraction: float,
            sample_group_column: str,
            control_group_name: str,
            experimental_group_name: str):

        self.normalized_intensity_table = normalized_intensity_table
        self.sample_table = sample_table
        self.probe_table = probe_table
        self.gene_sets_gmt = gene_sets_gmt
        self.gene_name_column = gene_name_column
        self.gene_description_column = gene_description_column
        self.heatmap_intensity_fraction = heatmap_intensity_fraction
        self.sample_group_column = sample_group_column
        self.control_group_name = control_group_name
        self.experimental_group_name = experimental_group_name

        self.read_tables()
        self.heatmap()
        self.pca()
        self.limma()
        self.gsea()
        self.clean_up()

    def read_tables(self):
        self.normalized_intensity_df = self.__read(self.normalized_intensity_table)
        self.sample_df = self.__read(self.sample_table)
        self.probe_df = self.__read(self.probe_table)

        for df in [self.normalized_intensity_df, self.sample_df, self.probe_df]:
            df.index.name = None  # make all final output files clean without index names

    def __read(self, file: str) -> pd.DataFrame:
        sep = ','
        for ext in ['.tsv', '.txt', '.tab']:
            if file.endswith(ext):
                sep = '\t'
                break
        return pd.read_csv(file, sep=sep, index_col=0)

    def heatmap(self):
        Heatmap(self.settings).main(
            normalized_intensity_df=self.normalized_intensity_df,
            heatmap_intensity_fraction=self.heatmap_intensity_fraction)

    def pca(self):
        PCA(self.settings).main(
            normalized_intensity_df=self.normalized_intensity_df,
            sample_df=self.sample_df,
            sample_group_column=self.sample_group_column)

    def limma(self):
        Limma(self.settings).main(
            normalized_intensity_df=self.normalized_intensity_df,
            sample_df=self.sample_df,
            sample_group_column=self.sample_group_column,
            control_group_name=self.control_group_name,
            experimental_group_name=self.experimental_group_name,
            probe_df=self.probe_df,
            gene_name_column=self.gene_name_column,
            gene_description_column=self.gene_description_column)

    def gsea(self):
        if self.gene_sets_gmt is None:
            return
        GSEA(self.settings).main(
            intensity_df=self.normalized_intensity_df,
            gene_info_df=self.probe_df,
            sample_info_df=self.sample_df,
            gene_name_column=self.gene_name_column,
            sample_group_column=self.sample_group_column,
            control_group_name=self.control_group_name,
            experimental_group_name=self.experimental_group_name,
            gene_sets_gmt=self.gene_sets_gmt)

    def clean_up(self):
        CleanUp(self.settings).main()


class CleanUp(Processor):

    def main(self):
        self.collect_files(file_ext='R', dstdir_name='r-scripts')
        self.collect_files(file_ext='log', dstdir_name='log')
        self.remove_workdir()

    def collect_files(self, file_ext: str, dstdir_name: str):
        files = get_files(
            source=self.outdir,
            endswith=file_ext)

        if len(files) == 0:
            return

        d = os.path.join(self.outdir, dstdir_name)
        os.makedirs(d, exist_ok=True)
        cmd = f'mv {self.outdir}/*.{file_ext} {d}/'
        self.call(cmd)

    def remove_workdir(self):
        if not self.debug:
            self.call(f'rm -r {self.workdir}')
