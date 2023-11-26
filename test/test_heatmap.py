import pandas as pd
from microarray_analysis.heatmap import Heatmap
from .setup import TestCase


class TestHeatmap(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        Heatmap(self.settings).main(
            normalized_intensity_df=pd.read_csv(f'{self.indir}/normalized_intensity.csv', index_col=0),
            probe_df=pd.read_csv(f'{self.indir}/probes.csv', index_col=0),
            gene_name_column='Gene Name',
            gene_description_column='Description',
            heatmap_intensity_fraction=0.8
        )
