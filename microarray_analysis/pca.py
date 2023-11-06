import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Tuple
from matplotlib.axes import Axes
from sklearn import decomposition
from .template import Processor


DSTDIR_NAME = 'pca'


class PCA(Processor):

    normalized_intensity_df: pd.DataFrame
    sample_df: pd.DataFrame
    sample_group_column: str

    def main(
            self,
            normalized_intensity_df: pd.DataFrame,
            sample_df: pd.DataFrame,
            sample_group_column: str):

        self.normalized_intensity_df = normalized_intensity_df
        self.sample_df = sample_df
        self.sample_group_column = sample_group_column

        OnePCA(self.settings).main(
            feature_by_sample_df=self.normalized_intensity_df,
            sample_df=self.sample_df,
            sample_group_column=self.sample_group_column,
            fname='pca'
        )


class OnePCA(Processor):

    XY_COLUMNS = ('PC 1', 'PC 2')

    feature_by_sample_df: pd.DataFrame  # row features x column samples
    sample_df: pd.DataFrame
    sample_group_column: str
    fname: str

    sample_coordinate_df: pd.DataFrame
    proportion_explained_series: pd.Series

    def main(
            self,
            feature_by_sample_df: pd.DataFrame,
            sample_df: pd.DataFrame,
            sample_group_column: str,
            fname: str):

        self.feature_by_sample_df = feature_by_sample_df
        self.sample_df = sample_df
        self.sample_group_column = sample_group_column
        self.fname = fname

        self.compute_pca()
        self.merge_group_info()
        self.make_dstdir()
        self.write_sample_coordinate()
        self.plot_sample_coordinate()
        self.write_proportion_explained()

    def compute_pca(self):
        self.sample_coordinate_df, self.proportion_explained_series = ComputePCA(self.settings).main(
            feature_by_sample_df=self.feature_by_sample_df)

    def merge_group_info(self):
        self.sample_coordinate_df = self.sample_coordinate_df.merge(
            right=self.sample_df,
            right_index=True,
            left_index=True,
            how='left'
        )

    def make_dstdir(self):
        os.makedirs(f'{self.outdir}/{DSTDIR_NAME}', exist_ok=True)

    def write_sample_coordinate(self):
        self.sample_coordinate_df.to_csv(
            f'{self.outdir}/{DSTDIR_NAME}/{self.fname}-sample-coordinate.csv'
        )

    def plot_sample_coordinate(self):
        ScatterPlot(self.settings).main(
            sample_coordinate_df=self.sample_coordinate_df,
            x_column=self.XY_COLUMNS[0],
            y_column=self.XY_COLUMNS[1],
            hue_column=self.sample_group_column,
            fname=f'{self.fname}-sample-coordinate'
        )

    def write_proportion_explained(self):
        self.proportion_explained_series.to_csv(
            f'{self.outdir}/{DSTDIR_NAME}/{self.fname}-proportion-explained.csv',
            header=['Proportion Explained']
        )


class ComputePCA(Processor):

    XY_COLUMNS = ('PC 1', 'PC 2')
    N_COMPONENTS = 2
    RANDOM_STATE = 1  # to ensure reproducible result

    feature_by_sample_df: pd.DataFrame

    sample_by_feature_df: pd.DataFrame
    embedding: decomposition.PCA
    sample_coordinate_df: pd.DataFrame
    proportion_explained_series: pd.Series

    def main(self, feature_by_sample_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.Series]:

        self.feature_by_sample_df = feature_by_sample_df

        self.transpose()
        self.set_embedding()
        self.fit_transform()
        self.set_proportion_explained_serise()

        return self.sample_coordinate_df, self.proportion_explained_series

    def transpose(self):
        self.sample_by_feature_df = self.feature_by_sample_df.transpose()

    def set_embedding(self):
        self.embedding = decomposition.PCA(
            n_components=self.N_COMPONENTS,
            copy=True,
            whiten=False,
            svd_solver='auto',
            tol=0.0,
            iterated_power='auto',
            random_state=self.RANDOM_STATE)

    def fit_transform(self):
        transformed = self.embedding.fit_transform(
            self.sample_by_feature_df.to_numpy()
        )

        self.sample_coordinate_df = pd.DataFrame(
            data=transformed,
            columns=self.XY_COLUMNS,
            index=self.sample_by_feature_df.index
        )

    def set_proportion_explained_serise(self):
        self.proportion_explained_series = pd.Series(self.embedding.explained_variance_ratio_)


class ScatterPlot(Processor):

    FIGSIZE = (8, 8)
    DPI = 600

    sample_coordinate_df: pd.DataFrame
    x_column: str
    y_column: str
    group_column: str
    fname: str

    ax: Axes

    def main(
            self,
            sample_coordinate_df: pd.DataFrame,
            x_column: str,
            y_column: str,
            hue_column: str,
            fname: str):

        self.sample_coordinate_df = sample_coordinate_df
        self.x_column = x_column
        self.y_column = y_column
        self.group_column = hue_column
        self.fname = fname

        self.init_figure()
        self.scatterplot()
        self.label_points()
        self.save_figure()

    def init_figure(self):
        plt.figure(figsize=self.FIGSIZE, dpi=self.DPI)

    def scatterplot(self):
        self.ax = sns.scatterplot(
            data=self.sample_coordinate_df,
            x=self.x_column,
            y=self.y_column,
            hue=self.group_column)

    def label_points(self):
        df = self.sample_coordinate_df
        for sample_name in df.index:
            self.ax.text(
                x=df.loc[sample_name, self.x_column],
                y=df.loc[sample_name, self.y_column],
                s=sample_name
            )

    def save_figure(self):
        for ext in ['pdf', 'png']:
            plt.savefig(f'{self.outdir}/{DSTDIR_NAME}/{self.fname}.{ext}', dpi=self.DPI)
        plt.close()
