import argparse
import microarray_analysis


__VERSION__ = '1.0.0'


PROG = 'python microarray_analysis'
DESCRIPTION = f'Affymetrix microarray analysis (version {__VERSION__}) by Yu-Cheng Lin (ylin@nycu.edu.tw)'
REQUIRED = [
    {
        'keys': ['-i', '--normalized-intensity-table'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the count table (gene rows x sample columns)',
        }
    },
    {
        'keys': ['-s', '--sample-table'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the sample info table (sample rows), the first (index) column should be sample ID',
        }
    },
    {
        'keys': ['-p', '--probe-table'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the gene info table (gene rows), the first (index) column should be gene ID',
        }
    },
]
OPTIONAL = [
    {
        'keys': ['-m', '--gene-sets-gmt'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'path to the gene sets gmt file for GSEA, if None then skip GSEA (default: %(default)s)',
        }
    },
    {
        'keys': ['--gene-name-column'],
        'properties': {
            'type': str,
            'required': False,
            'default': '"Gene Name"',
            'help': 'gene name (aka symbol) column in the gene-info-table (default: %(default)s)',
        }
    },
    {
        'keys': ['--gene-description-column'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'None',
            'help': 'gene description column in the gene-info-table, if None then no description will be used (default: %(default)s)',
        }
    },
    {
        'keys': ['--heatmap-intensity-fraction'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.8,
            'help': 'fraction of intensity to be included in the heatmap (default: %(default)s)',
        }
    },
    {
        'keys': ['--sample-group-column'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'Group',
            'help': 'sample group column in the sample-info-table (default: %(default)s)',
        }
    },
    {
        'keys': ['--control-group-name'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'WT',
            'help': 'control group name in the "sample group column" (default: %(default)s)',
        }
    },
    {
        'keys': ['--experimental-group-name'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'Mutant',
            'help': 'experimental group name in the "sample group column" (default: %(default)s)',
        }
    },
    {
        'keys': ['-o', '--outdir'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'rna_seq_analysis_outdir',
            'help': 'path to the output directory (default: %(default)s)',
        }
    },
    {
        'keys': ['-t', '--threads'],
        'properties': {
            'type': int,
            'required': False,
            'default': 4,
            'help': 'number of CPU threads (default: %(default)s)',
        }
    },
    {
        'keys': ['-d', '--debug'],
        'properties': {
            'action': 'store_true',
            'help': 'debug mode',
        }
    },
    {
        'keys': ['-h', '--help'],
        'properties': {
            'action': 'help',
            'help': 'show this help message',
        }
    },
    {
        'keys': ['-v', '--version'],
        'properties': {
            'action': 'version',
            'version': __VERSION__,
            'help': 'show version',
        }
    },
]


class EntryPoint:

    parser: argparse.ArgumentParser

    def main(self):
        self.set_parser()
        self.add_required_arguments()
        self.add_optional_arguments()
        self.run()

    def set_parser(self):
        self.parser = argparse.ArgumentParser(
            prog=PROG,
            description=DESCRIPTION,
            add_help=False,
            formatter_class=argparse.RawTextHelpFormatter)

    def add_required_arguments(self):
        group = self.parser.add_argument_group('required arguments')
        for item in REQUIRED:
            group.add_argument(*item['keys'], **item['properties'])

    def add_optional_arguments(self):
        group = self.parser.add_argument_group('optional arguments')
        for item in OPTIONAL:
            group.add_argument(*item['keys'], **item['properties'])

    def run(self):
        args = self.parser.parse_args()
        print(f'Start running RNA-seq Analysis version {__VERSION__}\n', flush=True)
        microarray_analysis.main(
            normalized_intensity_table=args.normalized_intensity_table,
            sample_table=args.sample_table,
            probe_table=args.probe_table,
            gene_sets_gmt=args.gene_sets_gmt,
            gene_name_column=args.gene_name_column,
            gene_description_column=args.gene_description_column,
            heatmap_intensity_fraction=args.heatmap_intensity_fraction,
            sample_group_column=args.sample_group_column,
            control_group_name=args.control_group_name,
            experimental_group_name=args.experimental_group_name,
            threads=args.threads,
            debug=args.debug,
            outdir=args.outdir)


if __name__ == '__main__':
    EntryPoint().main()
