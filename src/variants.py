#!/usr/bin/env python
"""
TETyper 2.0: Variant Calling Module
Handles SNP calling, VCF processing, and structural variant detection
"""

import os
import subprocess
import logging
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
import pysam
from collections import defaultdict, Counter


class ExternalToolWrapper:
    """Wrapper for external command execution with validation and error handling"""
    
    @staticmethod
    def validate_tool(tool_name: str) -> bool:
        """Check if external tool is available in PATH"""
        return shutil.which(tool_name) is not None
    
    @staticmethod
    def run_command(cmd: List[str], capture_output: bool = True, 
                   check_return: bool = True, cwd: Optional[str] = None) -> Tuple[int, str, str]:
        """Execute external command with proper error handling"""
        try:
            if capture_output:
                result = subprocess.run(cmd, capture_output=True, text=True, cwd=cwd)
            else:
                result = subprocess.run(cmd, cwd=cwd)
                
            if check_return and result.returncode != 0:
                raise RuntimeError(f"Command failed: {' '.join(cmd)}\nStderr: {result.stderr}")
                
            return result.returncode, getattr(result, 'stdout', ''), getattr(result, 'stderr', '')
            
        except Exception as e:
            raise RuntimeError(f"Failed to execute command {' '.join(cmd)}: {str(e)}")


class SNPCaller:
    """Handles SNP calling using samtools/bcftools pipeline"""
    
    def __init__(self, threads: int = 4, logger: Optional[logging.Logger] = None):
        self.threads = threads
        self.logger = logger or logging.getLogger(__name__)
        self.tool_wrapper = ExternalToolWrapper()
    
    def call_snps(self, bam_file: str, reference_fasta: str, output_vcf: str) -> str:
        """
        Call SNPs using samtools mpileup and bcftools call
        
        Args:
            bam_file: Input BAM file
            reference_fasta: Reference FASTA file
            output_vcf: Output VCF file path
            
        Returns:
            Path to output VCF file
        """
        for tool in ['samtools', 'bcftools']:
            if not self.tool_wrapper.validate_tool(tool):
                raise RuntimeError(f"{tool} not found in PATH")
        
        if not os.path.isfile(bam_file):
            raise FileNotFoundError(f"BAM file not found: {bam_file}")
        
        if not os.path.isfile(reference_fasta):
            raise FileNotFoundError(f"Reference file not found: {reference_fasta}")
        
        self.logger.info(f"Calling SNPs from {bam_file}")
        
        try:
            temp_bcf = f"{output_vcf}.tmp.bcf"
            
            cmd_mpileup = [
                'samtools', 'mpileup',
                '-u',
                '-A',
                '-I',
                '-f', reference_fasta,
                bam_file
            ]
            
            with open(temp_bcf, 'wb') as bcf_file:
                result_mpileup = subprocess.run(
                    cmd_mpileup, 
                    stdout=bcf_file, 
                    stderr=subprocess.PIPE,
                    text=False
                )
            
            if result_mpileup.returncode != 0:
                raise RuntimeError(f"samtools mpileup failed: {result_mpileup.stderr.decode()}")
            
            cmd_call = [
                'bcftools', 'call',
                '-m',
                '-v',
                '-o', output_vcf,
                temp_bcf
            ]
            
            returncode, stdout, stderr = self.tool_wrapper.run_command(cmd_call)
            
            if os.path.isfile(temp_bcf):
                os.remove(temp_bcf)
            
            if not os.path.isfile(output_vcf):
                raise RuntimeError("VCF file was not created")
            
            variant_count = self._count_variants_in_vcf(output_vcf)
            self.logger.info(f"SNP calling completed: {variant_count} variants in {output_vcf}")
            
            return output_vcf
            
        except Exception as e:
            for temp_file in [temp_bcf, output_vcf]:
                if os.path.isfile(temp_file):
                    os.remove(temp_file)
            raise e
    
    def _count_variants_in_vcf(self, vcf_file: str) -> int:
        """Count number of variants in VCF file"""
        count = 0
        try:
            with pysam.VariantFile(vcf_file, 'r') as vcf:
                for record in vcf:
                    count += 1
        except Exception:
            with open(vcf_file, 'r') as f:
                for line in f:
                    if not line.startswith('#') and line.strip():
                        count += 1
        return count


class VCFProcessor:
    """Processes VCF files to extract variant information"""
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        self.logger = logger or logging.getLogger(__name__)
    
    def parse_vcf(self, vcf_file: str, min_depth: int = 10, 
                 min_alt_depth: int = 3) -> Dict:
        """
        Parse VCF file and extract variant information
        
        Args:
            vcf_file: Input VCF file path
            min_depth: Minimum total depth for variant calling
            min_alt_depth: Minimum alternative allele depth
            
        Returns:
            Dictionary with variant information
        """
        if not os.path.isfile(vcf_file):
            raise FileNotFoundError(f"VCF file not found: {vcf_file}")
        
        variants = {
            'homozygous_snps': [],
            'heterozygous_snps': [],
            'homozygous_counts': [],
            'heterozygous_counts': [],
            'total_variants': 0,
            'filtered_variants': 0
        }
        
        try:
            with pysam.VariantFile(vcf_file, 'r') as vcf:
                for record in vcf:
                    variants['total_variants'] += 1
                    
                    if 'DP' in record.info:
                        total_depth = record.info['DP']
                    else:
                        total_depth = 0
                    
                    if total_depth < min_depth:
                        continue
                    
                    for sample in record.samples:
                        sample_data = record.samples[sample]
                        
                        if 'GT' not in sample_data or sample_data['GT'] is None:
                            continue
                        
                        gt = sample_data['GT']
                        
                        if 'AD' in sample_data and sample_data['AD'] is not None:
                            allele_depths = sample_data['AD']
                            ref_depth = allele_depths[0] if len(allele_depths) > 0 else 0
                            alt_depth = sum(allele_depths[1:]) if len(allele_depths) > 1 else 0
                        else:
                            ref_depth = 0
                            alt_depth = 0
                        
                        if alt_depth < min_alt_depth:
                            continue
                        
                        variants['filtered_variants'] += 1
                        
                        if gt == (1, 1):
                            snp_string = f"{record.ref}{record.pos}{record.alts[0]}"
                            variants['homozygous_snps'].append(snp_string)
                            variants['homozygous_counts'].append(f"{record.ref},{ref_depth},{record.alts[0]},{alt_depth}")
                            
                        elif gt == (0, 1) or gt == (1, 0):
                            ambiguous_base = self._get_ambiguous_base(record.ref, record.alts[0])
                            snp_string = f"{record.ref}{record.pos}{ambiguous_base}"
                            variants['heterozygous_snps'].append(snp_string)
                            variants['heterozygous_counts'].append(f"{record.ref},{ref_depth},{record.alts[0]},{alt_depth}")
        
        except Exception as e:
            self.logger.error(f"Error parsing VCF file {vcf_file}: {str(e)}")
            raise e
        
        self.logger.info(f"Parsed VCF: {variants['filtered_variants']} high-quality variants "
                        f"from {variants['total_variants']} total")
        
        return variants
    
    def _get_ambiguous_base(self, ref: str, alt: str) -> str:
        """Get IUPAC ambiguity code for heterozygous SNP"""
        bases = sorted([ref.upper(), alt.upper()])  
        ambiguity_map = {
            ('A', 'C'): 'M',
            ('A', 'G'): 'R', 
            ('A', 'T'): 'W',
            ('C', 'G'): 'S',
            ('C', 'T'): 'Y',
            ('G', 'T'): 'K'
        }
        
        key = tuple(bases)
        return ambiguity_map.get(key, 'N')


class StructuralVariantDetector:
    """Detects structural variants from assembly alignments"""
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        self.logger = logger or logging.getLogger(__name__)
        self.tool_wrapper = ExternalToolWrapper()
    
    def detect_deletions_from_blast(self, assembly_file: str, reference_fasta: str,
                                   blast_output: str) -> Dict:
        """
        Detect deletions by comparing assembly to reference via BLAST
        
        Args:
            assembly_file: Assembly FASTA file
            reference_fasta: Reference FASTA file
            blast_output: Output file for BLAST results
            
        Returns:
            Dictionary with deletion information
        """
        if not self.tool_wrapper.validate_tool('blastn'):
            raise RuntimeError("blastn not found in PATH")
        
        blast_db = f"{reference_fasta}.blastdb"
        self._create_blast_database(reference_fasta, blast_db)
        
        self._run_blast_alignment(assembly_file, blast_db, blast_output)
        
        deletions = self._parse_blast_for_deletions(blast_output, reference_fasta)
        
        return deletions
    
    def _create_blast_database(self, reference_fasta: str, db_name: str) -> None:
        """Create BLAST database from reference"""
        cmd = [
            'makeblastdb',
            '-dbtype', 'nucl',
            '-in', reference_fasta,
            '-out', db_name
        ]
        
        self.tool_wrapper.run_command(cmd)
        self.logger.debug(f"Created BLAST database: {db_name}")
    
    def _run_blast_alignment(self, query_file: str, blast_db: str, output_file: str) -> None:
        """Run BLAST alignment"""
        cmd = [
            'blastn',
            '-db', blast_db,
            '-query', query_file,
            '-outfmt', '6',
            '-out', output_file,
            '-num_threads', '4'
        ]
        
        self.tool_wrapper.run_command(cmd)
        self.logger.debug(f"BLAST alignment completed: {output_file}")
    
    def _parse_blast_for_deletions(self, blast_file: str, reference_fasta: str) -> Dict:
        """Parse BLAST results to identify deletions"""
        if not os.path.isfile(blast_file):
            return {'deletions': [], 'deletion_string': 'none'}
        
        from Bio import SeqIO
        ref_record = list(SeqIO.parse(reference_fasta, 'fasta'))[0]
        ref_length = len(ref_record.seq)
        
        hit_ranges = []
        with open(blast_file, 'r') as f:
            for line in f:
                if line.strip():
                    fields = line.strip().split('\t')
                    if len(fields) >= 10:
                        start = int(fields[8])
                        end = int(fields[9])  
                        hit_ranges.append((min(start, end), max(start, end)))
        
        if not hit_ranges:
            return {
                'deletions': [(1, ref_length)],
                'deletion_string': f"1-{ref_length}"
            }
        
        hit_ranges.sort()
        merged_ranges = []
        
        for start, end in hit_ranges:
            if not merged_ranges or start > merged_ranges[-1][1] + 1:
                merged_ranges.append((start, end))
            else:
                merged_ranges[-1] = (merged_ranges[-1][0], max(merged_ranges[-1][1], end))
        
        deletions = []
        last_end = 0
        
        for start, end in merged_ranges:
            if start > last_end + 1:
                deletions.append((last_end + 1, start - 1))
            last_end = end
        
        if last_end < ref_length:
            deletions.append((last_end + 1, ref_length))
        
        if deletions:
            deletion_string = ';'.join(f"{start}-{end}" for start, end in deletions)
        else:
            deletion_string = 'none'
        
        return {
            'deletions': deletions,
            'deletion_string': deletion_string,
            'covered_ranges': merged_ranges
        }


class VariantPipeline:
    """High-level variant calling pipeline"""
    
    def __init__(self, threads: int = 4, logger: Optional[logging.Logger] = None):
        self.threads = threads
        self.logger = logger or logging.getLogger(__name__)
        self.snp_caller = SNPCaller(threads, logger)
        self.vcf_processor = VCFProcessor(logger)
        self.sv_detector = StructuralVariantDetector(logger)
    
    def call_variants_complete(self, bam_file: str, reference_fasta: str,
                              assembly_file: Optional[str], output_dir: str,
                              sample_id: str, te_name: str) -> Dict:
        """
        Complete variant calling pipeline including SNPs and structural variants
        
        Args:
            bam_file: Input BAM file
            reference_fasta: Reference FASTA file
            assembly_file: Assembly file for structural variant detection
            output_dir: Output directory
            sample_id: Sample identifier
            te_name: TE name
            
        Returns:
            Dictionary with all variant information
        """
        os.makedirs(output_dir, exist_ok=True)
        
        results = {
            'sample_id': sample_id,
            'te_name': te_name,
            'snps': {},
            'structural_variants': {},
            'files': {}
        }
        
        vcf_file = os.path.join(output_dir, f"{sample_id}_{te_name}.vcf")
        
        try:
            self.snp_caller.call_snps(bam_file, reference_fasta, vcf_file)
            results['files']['vcf'] = vcf_file
            
            snp_data = self.vcf_processor.parse_vcf(vcf_file)
            results['snps'] = snp_data
            
        except Exception as e:
            self.logger.error(f"SNP calling failed for {sample_id} {te_name}: {str(e)}")
            results['snps']['error'] = str(e)
        
        if assembly_file and os.path.isfile(assembly_file):
            blast_file = os.path.join(output_dir, f"{sample_id}_{te_name}_blast.txt")
            
            try:
                sv_data = self.sv_detector.detect_deletions_from_blast(
                    assembly_file, reference_fasta, blast_file
                )
                results['structural_variants'] = sv_data
                results['files']['blast'] = blast_file
                
            except Exception as e:
                self.logger.error(f"Structural variant detection failed for {sample_id} {te_name}: {str(e)}")
                results['structural_variants']['error'] = str(e)
        
        return results


if __name__ == "__main__":
    import logging
    
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)
    
    wrapper = ExternalToolWrapper()
    
    required_tools = ['samtools', 'bcftools', 'blastn', 'makeblastdb']
    for tool in required_tools:
        if wrapper.validate_tool(tool):
            logger.info(f"{tool} is available")
        else:
            logger.error(f"{tool} is not available in PATH")
    
    snp_caller = SNPCaller(threads=4, logger=logger)
    vcf_processor = VCFProcessor(logger)
    sv_detector = StructuralVariantDetector(logger)
    pipeline = VariantPipeline(threads=4, logger=logger)
    
    logger.info("Variant calling module tests completed")