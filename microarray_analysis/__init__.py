import os
from .template import Settings
from .tools import get_temp_path
from .microarray_analysis import MicroarrayAnalysis


def main(
        normalized_intensity_table: str,
        sample_table: str,
        probe_table: str,
        gene_sets_gmt: str,
        gene_name_column: str,
        gene_description_column: str,
        heatmap_intensity_fraction: float,
        sample_group_column: str,
        control_group_name: str,
        experimental_group_name: str,
        threads: int,
        debug: bool,
        outdir: str):

    settings = Settings(
        workdir=get_temp_path(prefix='./microarray_analysis_workdir_'),
        outdir=outdir,
        threads=int(threads),
        debug=debug,
        mock=False)

    for d in [settings.workdir, settings.outdir]:
        os.makedirs(d, exist_ok=True)

    MicroarrayAnalysis(settings).main(
        normalized_intensity_table=normalized_intensity_table,
        sample_table=sample_table,
        probe_table=probe_table,
        gene_sets_gmt=None if gene_sets_gmt.lower() == 'none' else gene_sets_gmt,
        gene_name_column=gene_name_column,
        gene_description_column=None if gene_description_column.lower() == 'none' else gene_description_column,
        heatmap_intensity_fraction=heatmap_intensity_fraction,
        sample_group_column=sample_group_column,
        control_group_name=control_group_name,
        experimental_group_name=experimental_group_name)
