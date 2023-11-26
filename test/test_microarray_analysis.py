from .setup import TestCase
from microarray_analysis.microarray_analysis import MicroarrayAnalysis


class TestMicroarrayAnalysis(TestCase):

    def setUp(self):
        self.set_up(py_path=__file__)

    def tearDown(self):
        self.tear_down()

    def test_main(self):
        MicroarrayAnalysis(self.settings).main(
            normalized_intensity_table=f'{self.indir}/normalized_intensity.csv',
            sample_table=f'{self.indir}/samples.csv',
            probe_table=f'{self.indir}/probes.csv',
            gene_sets_gmt=f'{self.indir}/h.all.v2023.1.Hs.symbols.gmt',
            gene_name_column='Gene Name',
            gene_description_column='Description',
            heatmap_intensity_fraction=0.8,
            sample_group_column='Group',
            control_group_name='NaB',
            experimental_group_name='DEHP+NaB',
        )
