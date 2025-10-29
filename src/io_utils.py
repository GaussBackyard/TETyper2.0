#!/usr/bin/env python
"""
TETyper 2.0: Input/Output Module
Handles file validation, reference parsing, and output generation
"""

import os
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import json
import csv
import pandas as pd


class FileValidator:
    """Validates input files and directories"""
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        self.logger = logger or logging.getLogger(__name__)
    
    def validate_file_exists(self, filepath: str, description: str = "File") -> bool:
        """Check if a file exists and is readable"""
        if not os.path.isfile(filepath):
            raise FileNotFoundError(f"{description} not found: {filepath}")
        
        if os.path.getsize(filepath) == 0:
            raise ValueError(f"{description} is empty: {filepath}")
        
        return True
    
    def validate_fastq_pair(self, reads1: str, reads2: str) -> bool:
        """Validate paired-end FASTQ files"""
        self.validate_file_exists(reads1, "Forward reads file")
        self.validate_file_exists(reads2, "Reverse reads file")
        
        self.logger.info(f"Validated FASTQ pair: {reads1}, {reads2}")
        return True
    
    def validate_fasta_reference(self, fasta_path: str) -> Tuple[str, int]:
        """
        Validate FASTA reference and return contig name and length
        
        Returns:
            Tuple of (contig_name, sequence_length)
        """
        self.validate_file_exists(fasta_path, "FASTA reference")
        
        try:
            records = list(SeqIO.parse(fasta_path, "fasta"))
            if len(records) != 1:
                raise ValueError(f"Expected exactly one sequence in {fasta_path}, found {len(records)}")
            
            record = records[0]
            contig_name = record.id
            seq_length = len(record.seq)
            
            self.logger.info(f"Validated FASTA reference: {contig_name} ({seq_length} bp)")
            return contig_name, seq_length
            
        except Exception as e:
            raise ValueError(f"Error parsing FASTA file {fasta_path}: {str(e)}")


class ReferenceManager:
    """Manages reference sequences for different transposons"""
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        self.logger = logger or logging.getLogger(__name__)
        self.validator = FileValidator(logger)
    
    def create_ir_reference_files(self, te_name: str, irl_seq: str, irr_seq: str, 
                                 output_dir: str) -> Tuple[str, str]:
        """
        Create FASTA files for IRL and IRR sequences
        
        Args:
            te_name: Name of the transposon (e.g., 'Tn2', 'Tn4401')
            irl_seq: Inverted repeat left sequence
            irr_seq: Inverted repeat right sequence 
            output_dir: Directory to save FASTA files
            
        Returns:
            Tuple of (irl_fasta_path, irr_fasta_path)
        """
        os.makedirs(output_dir, exist_ok=True)
        
        irl_path = os.path.join(output_dir, f"{te_name}_IRL.fasta")
        with open(irl_path, 'w') as f:
            f.write(f">{te_name}_IRL\n{irl_seq}\n")
        
        irr_path = os.path.join(output_dir, f"{te_name}_IRR.fasta")
        with open(irr_path, 'w') as f:
            f.write(f">{te_name}_IRR\n{irr_seq}\n")
        
        self.validator.validate_fasta_reference(irl_path)
        self.validator.validate_fasta_reference(irr_path)
        
        self.logger.info(f"Created IR reference files for {te_name}")
        return irl_path, irr_path
    
    def load_reference_sequence(self, fasta_path: str) -> SeqRecord:
        """Load and return a reference sequence record"""
        self.validator.validate_fasta_reference(fasta_path)
        
        records = list(SeqIO.parse(fasta_path, "fasta"))
        return records[0]


class SampleSheetParser:
    """Parses sample sheets for batch processing"""
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        self.logger = logger or logging.getLogger(__name__)
        self.validator = FileValidator(logger)
    
    def parse_sample_sheet(self, sample_sheet_path: str) -> List[Dict[str, str]]:
        """
        Parse CSV/TSV sample sheet
        
        Expected format:
        sample_id,reads1,reads2
        CAV1026,/path/to/CAV1026_R1.fq.gz,/path/to/CAV1026_R2.fq.gz
        
        Returns:
            List of sample dictionaries
        """
        self.validator.validate_file_exists(sample_sheet_path, "Sample sheet")
        
        samples = []
        with open(sample_sheet_path, 'r') as f:
            dialect = csv.Sniffer().sniff(f.read(1024))
            f.seek(0)
            
            reader = csv.DictReader(f, delimiter=dialect.delimiter)
            
            required_columns = {'sample_id', 'reads1', 'reads2'}
            if not required_columns.issubset(set(reader.fieldnames)):
                raise ValueError(f"Sample sheet must contain columns: {required_columns}")
            
            for row_num, row in enumerate(reader, 1):
                try:
                    sample_id = row['sample_id'].strip()
                    reads1 = row['reads1'].strip()
                    reads2 = row['reads2'].strip()
                    
                    if not sample_id:
                        raise ValueError(f"Empty sample_id in row {row_num}")
                    
                    self.validator.validate_fastq_pair(reads1, reads2)
                    
                    samples.append({
                        'sample_id': sample_id,
                        'reads1': reads1,
                        'reads2': reads2
                    })
                    
                except Exception as e:
                    self.logger.warning(f"Skipping row {row_num} in sample sheet: {str(e)}")
                    continue
        
        if not samples:
            raise ValueError("No valid samples found in sample sheet")
        
        self.logger.info(f"Parsed {len(samples)} samples from sample sheet")
        return samples


class ResultsWriter:
    """Handles output generation in multiple formats"""
    
    def __init__(self, output_prefix: str, logger: Optional[logging.Logger] = None):
        self.output_prefix = output_prefix
        self.logger = logger or logging.getLogger(__name__)
    
    def write_json_results(self, results: Dict, suffix: str = "results") -> str:
        """Write results in JSON format"""
        output_path = f"{self.output_prefix}_{suffix}.json"
        
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        self.logger.info(f"JSON results written to {output_path}")
        return output_path
    
    def write_summary_tsv(self, results: Dict, suffix: str = "summary") -> str:
        """Write summary results in TSV format compatible with TETyper 1.0"""
        output_path = f"{self.output_prefix}_{suffix}.tsv"
        
        summary_rows = []
        
        if isinstance(results, dict) and 'sample_id' in results:
            summary_rows.append(self._extract_sample_summary(results))
        else:
            for sample_id, sample_result in results.items():
                if isinstance(sample_result, dict) and 'sample_id' in sample_result:
                    summary_rows.append(self._extract_sample_summary(sample_result))
        
        if summary_rows:
            df = pd.DataFrame(summary_rows)
            df.to_csv(output_path, sep='\t', index=False)
            self.logger.info(f"TSV summary written to {output_path}")
        
        return output_path
    
    def write_flanks_detail(self, results: Dict, suffix: str = "flanks") -> str:
        """Write detailed flanking sequence results"""
        output_path = f"{self.output_prefix}_{suffix}.tsv"
        
        flank_rows = []
        
        if isinstance(results, dict) and 'sample_id' in results:
            flank_rows.extend(self._extract_flank_details(results))
        else:
            for sample_id, sample_result in results.items():
                if isinstance(sample_result, dict):
                    flank_rows.extend(self._extract_flank_details(sample_result))
        
        if flank_rows:
            df = pd.DataFrame(flank_rows)
            df.to_csv(output_path, sep='\t', index=False)
            self.logger.info(f"Detailed flanks written to {output_path}")
        
        return output_path
    
    def _extract_sample_summary(self, sample_result: Dict) -> Dict:
        """Extract summary information for a single sample"""
        summary = {
            'Sample': sample_result.get('sample_id', 'Unknown'),
            'TE_Results': len(sample_result.get('te_results', {})),
        }
        
        te_results = sample_result.get('te_results', {})
        for te_name, te_result in te_results.items():
            left_flanks = len(te_result.get('left_flanks', []))
            right_flanks = len(te_result.get('right_flanks', []))
            
            summary[f'{te_name}_Left_Flanks'] = left_flanks
            summary[f'{te_name}_Right_Flanks'] = right_flanks
            summary[f'{te_name}_Total_Flanks'] = left_flanks + right_flanks
        
        return summary
    
    def _extract_flank_details(self, sample_result: Dict) -> List[Dict]:
        """Extract detailed flank information for a sample"""
        flank_details = []
        
        te_results = sample_result.get('te_results', {})
        for te_name, te_result in te_results.items():
            for flank_info in te_result.get('left_flanks', []):
                flank_details.append({
                    'Sample': sample_result.get('sample_id'),
                    'TE': te_name,
                    'Side': 'Left',
                    'Reference': flank_info.get('ref'),
                    'Sequence': flank_info.get('flank'),
                    'Count': flank_info.get('count'),
                    'Forward_Reads': flank_info.get('forward_reads'),
                    'Reverse_Reads': flank_info.get('reverse_reads'),
                    'Passes_Filter': flank_info.get('passes_filter')
                })
            
            for flank_info in te_result.get('right_flanks', []):
                flank_details.append({
                    'Sample': sample_result.get('sample_id'),
                    'TE': te_name,
                    'Side': 'Right',
                    'Reference': flank_info.get('ref'),
                    'Sequence': flank_info.get('flank'),
                    'Count': flank_info.get('count'),
                    'Forward_Reads': flank_info.get('forward_reads'),
                    'Reverse_Reads': flank_info.get('reverse_reads'),
                    'Passes_Filter': flank_info.get('passes_filter')
                })
        
        return flank_details
    
    def write_all_formats(self, results: Dict) -> Dict[str, str]:
        """Write results in all supported formats"""
        output_files = {}
        
        output_files['json'] = self.write_json_results(results)
        
        output_files['summary'] = self.write_summary_tsv(results)
        
        output_files['flanks'] = self.write_flanks_detail(results)
        
        return output_files


class DirectoryManager:
    """Manages working directories and temporary files"""
    
    def __init__(self, base_output_dir: str, logger: Optional[logging.Logger] = None):
        self.base_output_dir = Path(base_output_dir)
        self.logger = logger or logging.getLogger(__name__)
    
    def create_working_directory(self, sample_id: str) -> str:
        """Create and return a working directory for a sample"""
        sample_dir = self.base_output_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        
        self.logger.debug(f"Created working directory: {sample_dir}")
        return str(sample_dir)
    
    def get_temp_file_path(self, sample_id: str, filename: str) -> str:
        """Get path for temporary file in sample directory"""
        sample_dir = self.base_output_dir / sample_id
        return str(sample_dir / filename)
    
    def cleanup_temp_files(self, sample_id: str, keep_files: Optional[List[str]] = None) -> None:
        """Clean up temporary files for a sample"""
        keep_files = keep_files or []
        sample_dir = self.base_output_dir / sample_id
        
        if sample_dir.exists():
            for file_path in sample_dir.iterdir():
                if file_path.name not in keep_files:
                    try:
                        file_path.unlink()
                        self.logger.debug(f"Removed temporary file: {file_path}")
                    except Exception as e:
                        self.logger.warning(f"Failed to remove {file_path}: {str(e)}")


if __name__ == "__main__":
    import logging
    
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    
    validator = FileValidator(logger)
    
    ref_manager = ReferenceManager(logger)
    
    try:
        irl_path, irr_path = ref_manager.create_ir_reference_files(
            te_name="Tn2",
            irl_seq="GGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGG",
            irr_seq="GGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGCA", 
            output_dir="/home/jl9gx/plasmid/TETyper2.0/tetyper2_refs"
        )
        logger.info(f"Created reference files: {irl_path}, {irr_path}")
        
    except Exception as e:
        logger.error(f"Failed to create reference files: {str(e)}")
    
    writer = ResultsWriter("/home/jl9gx/plasmid/TETyper2.0/test_output")
    
    example_results = {
        "CAV1026": {
            "sample_id": "CAV1026",
            "te_results": {
                "Tn2": {
                    "left_flanks": [
                        {"ref": "Tn2_IRL", "flank": "ATCGATCGAT", "count": 15, 
                         "forward_reads": 8, "reverse_reads": 7, "passes_filter": True}
                    ],
                    "right_flanks": [
                        {"ref": "Tn2_IRR", "flank": "GCTAGCTAGT", "count": 12,
                         "forward_reads": 6, "reverse_reads": 6, "passes_filter": True}
                    ]
                }
            }
        }
    }
    
    try:
        output_files = writer.write_all_formats(example_results)
        logger.info(f"Generated output files: {output_files}")
    except Exception as e:
        logger.error(f"Failed to write results: {str(e)}")