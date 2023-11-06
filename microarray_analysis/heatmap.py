import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Tuple
from .template import Processor


class Heatmap(Processor):

    normalized_intensity_df: pd.DataFrame
    heatmap_intensity_fraction: float

    def main(
            self,
            normalized_intensity_df: pd.DataFrame,
            heatmap_intensity_fraction: float):

        self.normalized_intensity_df = normalized_intensity_df
        self.heatmap_intensity_fraction = heatmap_intensity_fraction

        OneHeatmap(self.settings).main(
            feature_by_sample_df=self.normalized_intensity_df,
            heatmap_intensity_fraction=self.heatmap_intensity_fraction,
            fname='heatmap'
        )


class OneHeatmap(Processor):

    df: pd.DataFrame
    heatmap_intensity_fraction: float
    fname: str

    def main(
            self,
            feature_by_sample_df: pd.DataFrame,
            heatmap_intensity_fraction: float,
            fname: str):

        self.df = feature_by_sample_df.copy()
        self.heatmap_intensity_fraction = heatmap_intensity_fraction
        self.fname = fname

        self.filter_by_cumulative_intensity()
        self.df = np.log10(self.df)
        self.clustermap()

    def filter_by_cumulative_intensity(self):
        self.logger.info(f'Filter by cumulative intensity for "{self.fname}"')
        self.df = FilterByCumulativeReads(self.settings).main(
            df=self.df,
            heatmap_intensity_fraction=self.heatmap_intensity_fraction)

    def clustermap(self):
        Clustermap(self.settings).main(
            data=self.df,
            fname=self.fname)


class FilterByCumulativeReads(Processor):

    ROW_SUM = 'Row Sum'
    CUMULATIVE_SUM = 'Cumulative Sum'

    df: pd.DataFrame
    heatmap_read_fraction: float

    def main(
            self,
            df: pd.DataFrame,
            heatmap_intensity_fraction: float) -> pd.DataFrame:

        self.df = df
        self.heatmap_read_fraction = heatmap_intensity_fraction

        assert 0 < self.heatmap_read_fraction < 1

        self.sum_each_row()
        self.sort_by_sum()
        self.cumulate_sum()
        self.divide_by_total()
        self.filter()
        self.clean_up()

        return self.df

    def sum_each_row(self):
        self.df[self.ROW_SUM] = np.sum(self.df.to_numpy(), axis=1)

    def sort_by_sum(self):
        self.df = self.df.sort_values(
            by=self.ROW_SUM,
            ascending=False)

    def cumulate_sum(self):
        l = []
        c = 0
        for s in self.df[self.ROW_SUM]:
            c += s
            l.append(c)
        self.df[self.CUMULATIVE_SUM] = l

    def divide_by_total(self):
        total = self.df[self.CUMULATIVE_SUM][-1]
        self.df[self.CUMULATIVE_SUM] = self.df[self.CUMULATIVE_SUM] / total

    def filter(self):
        total_rows = len(self.df)
        n_rows = total_rows
        for i, v in enumerate(self.df[self.CUMULATIVE_SUM]):
            if v > self.heatmap_read_fraction:
                n_rows = i
                break
        self.df = self.df.iloc[:n_rows, ]

        msg = f'''\
Total genes: {total_rows}
For heatmap, keep the most abundant genes covering {self.heatmap_read_fraction * 100:.2f}% of counts: {n_rows} genes'''
        self.logger.info(msg)

    def clean_up(self):
        self.df = self.df.drop(
            columns=[self.ROW_SUM, self.CUMULATIVE_SUM]
        )


class Clustermap(Processor):

    CLUSTER_COLUMNS = True
    COLORMAP = 'PuBu'
    Y_LABEL_CHAR_WIDTH = 0.08
    X_LABEL_CHAR_WIDTH = 0.08
    XTICKLABELS = True
    YTICKLABELS = False
    CELL_WIDTH = 0.3
    CELL_HEIGHT = 0.01
    LINEWIDTH = 0.
    DENDROGRAM_RATIO = (0.1, 0.05)
    COLORBAR_WIDTH = 0.01
    COLORBAR_HORIZONTAL_POSITION = 1.
    DPI = 600
    DSTDIR_NAME = 'heatmap'

    data: pd.DataFrame
    fname: str

    x_label_padding: float
    y_label_padding: float
    figsize: Tuple[float, float]
    grid: sns.matrix.ClusterGrid

    def main(self, data: pd.DataFrame, fname: str):
        self.data = data
        self.fname = fname

        self.set_figsize()
        self.clustermap()
        self.config_clustermap()
        self.make_dstdir()
        self.save_fig()
        self.save_csv()

    def set_figsize(self):
        self.__set_x_y_label_padding()
        w = (len(self.data.columns) * self.CELL_WIDTH) + self.y_label_padding
        h = (len(self.data.index) * self.CELL_HEIGHT) + self.x_label_padding
        self.figsize = (w, h)

    def __set_x_y_label_padding(self):
        max_x_label_length = pd.Series(self.data.columns).apply(len).max()
        self.x_label_padding = max_x_label_length * self.X_LABEL_CHAR_WIDTH

        max_y_label_length = pd.Series(self.data.index).apply(len).max()
        self.y_label_padding = max_y_label_length * self.Y_LABEL_CHAR_WIDTH

    def clustermap(self):
        self.grid = sns.clustermap(
            data=self.data,
            cmap=self.COLORMAP,
            figsize=self.figsize,
            xticklabels=self.XTICKLABELS,
            yticklabels=self.YTICKLABELS,
            col_cluster=self.CLUSTER_COLUMNS,
            dendrogram_ratio=self.DENDROGRAM_RATIO,
            linewidth=self.LINEWIDTH)
        self.__set_plotted_data()

    def __set_plotted_data(self):
        self.data = self.grid.__dict__['data2d']

    def config_clustermap(self):
        self.__set_x_y_axes()
        self.__set_colorbar()

    def __set_x_y_axes(self):
        heatmap = self.grid.ax_heatmap

        n_rows = len(self.data)
        heatmap.set_ylim(n_rows, 0)  # first and last row not be chopped half

        heatmap.tick_params(
            axis='x',
            bottom=True,
            top=False,
            labelbottom=True,
            labeltop=False,
            labelrotation=90
        )
        heatmap.tick_params(
            axis='y',
            left=False,
            right=True,
            labelleft=False,
            labelright=True,
            labelrotation=0
        )

    def __set_colorbar(self):
        colorbar = self.grid.cax
        p = colorbar.get_position()
        colorbar.set_position([
            self.COLORBAR_HORIZONTAL_POSITION,  # left
            p.y0,  # bottom
            self.COLORBAR_WIDTH,  # width
            p.height  # height
        ])

    def make_dstdir(self):
        os.makedirs(f'{self.outdir}/{self.DSTDIR_NAME}', exist_ok=True)

    def save_fig(self):
        # must use grid.savefig(), but not plt.savefig()
        # plt.savefig() crops out the colorbar

        dpi = self.__downsize_dpi_if_too_large()
        for ext in ['pdf', 'png']:
            self.grid.savefig(f'{self.outdir}/{self.DSTDIR_NAME}/{self.fname}.{ext}', dpi=dpi)
        plt.close()

    def __downsize_dpi_if_too_large(self) -> int:
        longer_side = max(self.figsize)
        dpi = self.DPI
        while longer_side * dpi >= 2**16:  # 2^16 is the limit
            dpi = int(dpi/2)  # downsize
        return dpi

    def save_csv(self):
        self.data.to_csv(f'{self.outdir}/{self.DSTDIR_NAME}/{self.fname}.csv', index=True)
