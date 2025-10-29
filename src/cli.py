#!/usr/bin/env python
"""
TETyper 2.0: Command Line Interface Module
Handles argument parsing and main pipeline coordination
"""

import argparse
import sys
import os
import logging
from pathlib import Path
from typing import Dict, List, Optional
import json

from io_utils import FileValidator, ReferenceManager, SampleSheetParser, ResultsWriter, DirectoryManager
from mapping import TEMapper, ParallelMapper
from assembly import AssemblyPipeline
from variants import VariantPipeline
from flank_extractor import FlankExtractor, ProcessingParameters

VERSION = "2.0.0-dev"
DESCRIPTION = "TETyper 2.0: Enhanced transposable element typing with multi-TE support"

DEFAULT_IR_SEQUENCES = {
    "Tn2": {
        "IRL": "GGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGG",
        "IRR": "GGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGCA"
    },
    "Tn4401": {
        "IRL": "GGGGTTCTAAGCGGGAATCCCAGAAAATTCCGTCATTCCG", 
        "IRR": "GGGGGGGTAAGCGGGAECCCCAGAAAATTCCGCCATTCCG"
    }
}


class TETyper2Pipeline:
    """Main TETyper 2.0 pipeline coordinator"""
    
    def __init__(self, args):
        self.args = args
        self.logger = self._setup_logging()
        
        self.validator = FileValidator(self.logger)
        self.ref_manager = ReferenceManager(self.logger)
        self.dir_manager = DirectoryManager(args.output_dir)
        self.results_writer = ResultsWriter(args.output_prefix, self.logger)
        
        self.processing_params = ProcessingParameters(
            flank_length=args.flank_length,
            min_mapped_length=args.min_mapped_length,
            min_quality=args.min_quality,
            min_reads=args.min_reads,
            min_each_strand=args.min_each_strand,
            threads=args.threads
        )
        
        self.mapper = TEMapper(threads=args.threads, logger=self.logger)
        self.parallel_mapper = ParallelMapper(
            max_workers=args.max_workers,
            threads_per_worker=args.threads,
            logger=self.logger
        )
        self.assembly_pipeline = AssemblyPipeline(
            threads=args.threads,
            memory_limit=args.memory_limit,
            logger=self.logger
        )
        self.variant_pipeline = VariantPipeline(threads=args.threads, logger=self.logger)
        self.flank_extractor = FlankExtractor(self.processing_params, self.logger)
        
        self._validate_dependencies()
    
    def _setup_logging(self) -> logging.Logger:
        """Setup logging configuration"""
        log_level = logging.DEBUG if self.args.verbose else logging.INFO
        logging.basicConfig(
            level=log_level,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.StreamHandler(),
                logging.FileHandler(f"{self.args.output_prefix}.log")
            ]
        )
        
        logger = logging.getLogger('TETyper2')
        logger.info(f"TETyper 2.0 v{VERSION} started")
        logger.info(f"Command: {' '.join(sys.argv)}")
        
        return logger
    
    def _validate_dependencies(self):
        """Validate that all required external tools are available"""
        required_tools = ['bwa', 'samtools', 'bcftools']
        
        if self.args.enable_assembly:
            required_tools.append('spades.py')
        
        if self.args.enable_variants:
            required_tools.extend(['blastn', 'makeblastdb'])
        
        missing_tools = []
        from mapping import ExternalToolWrapper
        
        for tool in required_tools:
            if not ExternalToolWrapper.validate_tool(tool):
                missing_tools.append(tool)
        
        if missing_tools:
            raise RuntimeError(f"Missing required tools: {', '.join(missing_tools)}")
        
        self.logger.info("All required external tools are available")
    
    def load_te_configurations(self) -> List[Dict]:
        """Load TE configurations from arguments"""
        te_configs = []
        
        if self.args.te_config:
            for config_file in self.args.te_config:
                self.validator.validate_file_exists(config_file, "TE configuration file")
                with open(config_file, 'r') as f:
                    config_data = json.load(f)
                    te_configs.append(config_data)
        
        if self.args.te_name:
            for te_name in self.args.te_name:
                if te_name in DEFAULT_IR_SEQUENCES:
                    ir_seqs = DEFAULT_IR_SEQUENCES[te_name]
                    te_configs.append({
                        'name': te_name,
                        'irl_sequence': ir_seqs['IRL'],
                        'irr_sequence': ir_seqs['IRR']
                    })
                else:
                    raise ValueError(f"Unknown built-in TE: {te_name}")
        
        if not te_configs:
            raise ValueError("No TE configurations specified")
        
        self.logger.info(f"Loaded {len(te_configs)} TE configurations")
        return te_configs
    
    def prepare_references(self, te_configs: List[Dict]) -> Dict:
        """Prepare reference files for all TEs"""
        ref_files = {}
        
        for te_config in te_configs:
            te_name = te_config['name']
            
            if 'reference_fasta' in te_config and te_config['reference_fasta']:
                ref_files[te_name] = {
                    'reference': te_config['reference_fasta']
                }
            elif 'irl_sequence' in te_config and 'irr_sequence' in te_config:
                refs_dir = os.path.join(self.args.output_dir, 'references')
                
                irl_path, irr_path = self.ref_manager.create_ir_reference_files(
                    te_name=te_name,
                    irl_seq=te_config['irl_sequence'],
                    irr_seq=te_config['irr_sequence'],
                    output_dir=refs_dir
                )
                
                ref_files[te_name] = {
                    'irl_fasta': irl_path,
                    'irr_fasta': irr_path
                }
            else:
                raise ValueError(f"TE {te_name} requires either reference_fasta or IRL/IRR sequences")
        
        return ref_files
    
    def run_single_sample(self, sample_info: Dict, te_configs: List[Dict], 
                         ref_files: Dict) -> Dict:
        """Process a single sample through the complete pipeline"""
        sample_id = sample_info['sample_id']
        reads1 = sample_info['reads1']
        reads2 = sample_info['reads2']
        
        self.logger.info(f"Processing sample: {sample_id}")
        
        sample_dir = self.dir_manager.create_working_directory(sample_id)
        
        results = {
            'sample_id': sample_id,
            'te_results': {},
            'summary': {}
        }
        
        for te_config in te_configs:
            te_name = te_config['name']
            
            try:
                self.logger.info(f"Processing {sample_id} for {te_name}")
                
                te_result = self._process_te_for_sample(
                    sample_id, reads1, reads2, te_name,
                    ref_files[te_name], sample_dir
                )
                
                results['te_results'][te_name] = te_result
                
            except Exception as e:
                self.logger.error(f"Failed to process {te_name} for {sample_id}: {str(e)}")
                results['te_results'][te_name] = {'error': str(e)}
        
        results['summary'] = self._generate_sample_summary(results['te_results'])
        
        return results
    
    def _process_te_for_sample(self, sample_id: str, reads1: str, reads2: str,
                              te_name: str, ref_files: Dict, sample_dir: str) -> Dict:
        """Process a single TE for a sample"""
        
        te_result = {
            'te_name': te_name,
            'left_flanks': [],
            'right_flanks': [],
            'variants': {},
            'assembly_stats': {}
        }
        
        if 'irl_fasta' in ref_files and 'irr_fasta' in ref_files:
            irl_bam, irr_bam = self.mapper.map_sample_to_te(
                sample_id=sample_id,
                reads1=reads1,
                reads2=reads2,
                irl_fasta=ref_files['irl_fasta'],
                irr_fasta=ref_files['irr_fasta'],
                output_dir=sample_dir
            )
            
            flank_results = self.flank_extractor.extract_flanks_for_te(
                sample_id=sample_id,
                bam_left=irl_bam,
                bam_right=irr_bam,
                irl_fasta=ref_files['irl_fasta'],
                irr_fasta=ref_files['irr_fasta'],
                te_name=te_name
            )
            
            te_result['left_flanks'] = flank_results.get('left_flanks', [])
            te_result['right_flanks'] = flank_results.get('right_flanks', [])
            
            primary_bam = irl_bam
            reference_fasta = ref_files['irl_fasta']
            
        else:
            reference_fasta = ref_files['reference']
            primary_bam = self.mapper.map_reads_to_ir(
                sample_id=sample_id,
                reads1=reads1,
                reads2=reads2,
                ir_fasta=reference_fasta,
                output_dir=sample_dir,
                ir_type="reference"
            )
        
        assembly_file = None
        if self.args.enable_assembly:
            try:
                assembly_file, assembly_stats = self.assembly_pipeline.assemble_from_bam(
                    bam_file=primary_bam,
                    output_dir=sample_dir,
                    sample_id=sample_id,
                    te_name=te_name,
                    keep_intermediate=self.args.keep_intermediate
                )
                te_result['assembly_stats'] = assembly_stats
            except Exception as e:
                self.logger.warning(f"Assembly failed for {sample_id} {te_name}: {str(e)}")
        
        if self.args.enable_variants:
            try:
                variant_results = self.variant_pipeline.call_variants_complete(
                    bam_file=primary_bam,
                    reference_fasta=reference_fasta,
                    assembly_file=assembly_file,
                    output_dir=sample_dir,
                    sample_id=sample_id,
                    te_name=te_name
                )
                te_result['variants'] = variant_results
            except Exception as e:
                self.logger.warning(f"Variant calling failed for {sample_id} {te_name}: {str(e)}")
        
        return te_result
    
    def _generate_sample_summary(self, te_results: Dict) -> Dict:
        """Generate summary statistics for a sample"""
        summary = {
            'total_tes_processed': len(te_results),
            'successful_tes': sum(1 for result in te_results.values() if 'error' not in result),
            'total_left_flanks': 0,
            'total_right_flanks': 0
        }
        
        for te_name, te_result in te_results.items():
            if 'error' not in te_result:
                summary['total_left_flanks'] += len(te_result.get('left_flanks', []))
                summary['total_right_flanks'] += len(te_result.get('right_flanks', []))
        
        return summary
    
    def run_batch(self, sample_sheet_path: str, te_configs: List[Dict], 
                 ref_files: Dict) -> Dict:
        """Process multiple samples in batch mode"""
        parser = SampleSheetParser(self.logger)
        samples = parser.parse_sample_sheet(sample_sheet_path)
        
        self.logger.info(f"Processing {len(samples)} samples in batch mode")
        
        results = {}
        
        if self.args.max_workers > 1:
            self.logger.info(f"Using parallel processing with {self.args.max_workers} workers")
            
            for sample in samples:
                try:
                    result = self.run_single_sample(sample, te_configs, ref_files)
                    results[sample['sample_id']] = result
                except Exception as e:
                    self.logger.error(f"Failed to process sample {sample['sample_id']}: {str(e)}")
                    results[sample['sample_id']] = {'error': str(e)}
        else:
            for sample in samples:
                try:
                    result = self.run_single_sample(sample, te_configs, ref_files)
                    results[sample['sample_id']] = result
                except Exception as e:
                    self.logger.error(f"Failed to process sample {sample['sample_id']}: {str(e)}")
                    results[sample['sample_id']] = {'error': str(e)}
        
        return results
    
    def run(self):
        """Main execution method"""
        try:
            te_configs = self.load_te_configurations()
            
            ref_files = self.prepare_references(te_configs)
            
            if self.args.sample_sheet:
                results = self.run_batch(self.args.sample_sheet, te_configs, ref_files)
            else:
                sample_info = {
                    'sample_id': self.args.sample_id or Path(self.args.reads1).stem,
                    'reads1': self.args.reads1,
                    'reads2': self.args.reads2
                }
                results = self.run_single_sample(sample_info, te_configs, ref_files)
            
            if self.args.output_format in ['json', 'both']:
                self.results_writer.write_json_results(results)
            
            if self.args.output_format in ['tsv', 'both']:
                self.results_writer.write_summary_tsv(results)
                self.results_writer.write_flanks_detail(results)
            
            self.logger.info("TETyper 2.0 completed successfully")
            return 0
            
        except Exception as e:
            self.logger.error(f"TETyper 2.0 failed: {str(e)}")
            return 1


def setup_argument_parser() -> argparse.ArgumentParser:
    """Setup command line argument parser"""
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--version', action='version', version=f'TETyper {VERSION}')
    
    input_group = parser.add_argument_group('Input options')
    input_group.add_argument('-1', '--reads1', type=str,
                           help='Forward reads FASTQ file')
    input_group.add_argument('-2', '--reads2', type=str,
                           help='Reverse reads FASTQ file')
    input_group.add_argument('--sample-sheet', type=str,
                           help='CSV/TSV file with sample information for batch processing')
    input_group.add_argument('--sample-id', type=str,
                           help='Sample identifier (default: derived from reads1 filename)')
    
    te_group = parser.add_argument_group('Transposable element configuration')
    te_group.add_argument('--te-config', type=str, action='append',
                         help='JSON configuration file for TE (can be specified multiple times)')
    te_group.add_argument('--te-name', type=str, action='append',
                         choices=list(DEFAULT_IR_SEQUENCES.keys()),
                         help='TE name from built-in configurations')
    
    proc_group = parser.add_argument_group('Processing parameters')
    proc_group.add_argument('--flank-length', type=int, default=10,
                          help='Length of flanking sequences to extract')
    proc_group.add_argument('--min-mapped-length', type=int, default=5,
                          help='Minimum mapped length for flank extraction')
    proc_group.add_argument('--min-reads', type=int, default=10,
                          help='Minimum read support for flanking sequences')
    proc_group.add_argument('--min-each-strand', type=int, default=3,
                          help='Minimum reads per strand for flanking sequences')
    proc_group.add_argument('--min-quality', type=int, default=20,
                          help='Minimum base quality for flanking sequences')
    proc_group.add_argument('--threads', type=int, default=4,
                          help='Number of threads for external tools')
    proc_group.add_argument('--max-workers', type=int, default=1,
                          help='Number of parallel workers for batch processing')
    proc_group.add_argument('--memory-limit', type=int, default=8,
                          help='Memory limit in GB for SPAdes assembly')
    
    output_group = parser.add_argument_group('Output options')
    output_group.add_argument('-o', '--output-prefix', type=str, required=True,
                            help='Output prefix for all generated files')
    output_group.add_argument('--output-dir', type=str, default='.',
                            help='Output directory for intermediate and final files')
    output_group.add_argument('--output-format', type=str, default='both',
                            choices=['json', 'tsv', 'both'],
                            help='Output format')
    
    pipeline_group = parser.add_argument_group('Pipeline options')
    pipeline_group.add_argument('--enable-assembly', action='store_true',
                               help='Enable assembly-based structural variant detection')
    pipeline_group.add_argument('--enable-variants', action='store_true',
                               help='Enable SNP and structural variant calling')
    pipeline_group.add_argument('--keep-intermediate', action='store_true',
                               help='Keep intermediate files (BAM, FASTQ, etc.)')
    
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose logging')
    
    return parser


def validate_arguments(args):
    """Validate command line arguments"""
    if not args.sample_sheet:
        if not args.reads1 or not args.reads2:
            raise ValueError("Either --sample-sheet or both --reads1 and --reads2 must be provided")
    
    if not args.te_config and not args.te_name:
        raise ValueError("At least one TE configuration must be specified via --te-config or --te-name")
    
    validator = FileValidator()
    
    if args.reads1 and args.reads2:
        validator.validate_fastq_pair(args.reads1, args.reads2)
    
    if args.sample_sheet:
        validator.validate_file_exists(args.sample_sheet, "Sample sheet")
    
    os.makedirs(args.output_dir, exist_ok=True)


def main():
    """Main entry point"""
    parser = setup_argument_parser()
    args = parser.parse_args()
    
    try:
        validate_arguments(args)
        
        pipeline = TETyper2Pipeline(args)
        exit_code = pipeline.run()
        
        sys.exit(exit_code)
        
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()