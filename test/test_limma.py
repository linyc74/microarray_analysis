import pandas as pd
from microarray_analysis.limma import Limma
from .setup import TestCase


class TestLimma(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        Limma(self.settings).main(
            normalized_intensity_df=pd.read_csv(f'{self.indir}/normalized_intensity.csv', index_col=0),
            sample_df=pd.read_csv(f'{self.indir}/samples.csv', index_col=0),
            sample_group_column='Group',
            control_group_name='NaB',
            experimental_group_name='DEHP+NaB',
            probe_df=pd.read_csv(f'{self.indir}/probes.csv', index_col=0),
            gene_name_column='Gene Name',
            gene_description_column='Description',
        )
