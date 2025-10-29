#!/usr/bin/env python
"""
TETyper 2.0: Mapping Module
Handles BWA alignment and BAM file processing with consistent contig naming
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed
import shutil


class ExternalToolWrapper:
    """Wrapper for external command execution with validation and error handling"""
    
    @staticmethod
    def validate_tool(tool_name: str) -> bool:
        """Check if external tool is available in PATH"""
        return shutil.which(tool_name) is not None
    
    @staticmethod
    def run_command(cmd: List[str], capture_output: bool = True, 
                   check_return: bool = True, cwd: Optional[str] = None) -> Tuple[int, str, str]:
        """
        Execute external command with proper error handling
        
        Args:
            cmd: Command list
            capture_output: Whether to capture stdout/stderr
            check_return: Whether to check return code
            cwd: Working directory for command
            
        Returns:
            Tuple of (return_code, stdout, stderr)
            
        Raises:
            RuntimeError: If command fails and check_return is True
        """
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


class BWAIndexer:
    """Handles BWA index creation for reference sequences"""
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        self.logger = logger or logging.getLogger(__name__)
        self.tool_wrapper = ExternalToolWrapper()
    
    def create_index(self, reference_fasta: str) -> bool:
        """
        Create BWA index for reference FASTA file
        
        Args:
            reference_fasta: Path to FASTA reference file
            
        Returns:
            True if indexing successful
        """
        if not os.path.isfile(reference_fasta):
            raise FileNotFoundError(f"Reference file not found: {reference_fasta}")
        
        if not self.tool_wrapper.validate_tool('bwa'):
            raise RuntimeError("BWA not found in PATH")
        
        index_files = [f"{reference_fasta}.{ext}" for ext in ['amb', 'ann', 'bwt', 'pac', 'sa']]
        if all(os.path.isfile(idx_file) for idx_file in index_files):
            self.logger.info(f"BWA index already exists for {reference_fasta}")
            return True
        
        self.logger.info(f"Creating BWA index for {reference_fasta}")
        
        cmd = ['bwa', 'index', reference_fasta]
        returncode, stdout, stderr = self.tool_wrapper.run_command(cmd)
        
        if all(os.path.isfile(idx_file) for idx_file in index_files):
            self.logger.info(f"BWA index created successfully for {reference_fasta}")
            return True
        else:
            raise RuntimeError(f"BWA index creation failed for {reference_fasta}")


class BWAMapper:
    """Handles BWA alignment operations"""
    
    def __init__(self, threads: int = 4, logger: Optional[logging.Logger] = None):
        self.threads = threads
        self.logger = logger or logging.getLogger(__name__)
        self.tool_wrapper = ExternalToolWrapper()
        self.indexer = BWAIndexer(logger)
    
    def align_reads(self, reference_fasta: str, reads1: str, reads2: str, 
                   output_sam: str, sample_id: str = None) -> str:
        """
        Align paired-end reads to reference using BWA MEM
        
        Args:
            reference_fasta: Path to reference FASTA file
            reads1: Path to forward reads FASTQ
            reads2: Path to reverse reads FASTQ  
            output_sam: Path for output SAM file
            sample_id: Sample ID for read group
            
        Returns:
            Path to output SAM file
        """
        self.indexer.create_index(reference_fasta)
        
        cmd = [
            'bwa', 'mem',
            '-t', str(self.threads),
        ]
        
        if sample_id:
            rg_string = f"@RG\\tID:{sample_id}\\tSM:{sample_id}\\tPL:ILLUMINA"
            cmd.extend(['-R', rg_string])
        
        cmd.extend([reference_fasta, reads1, reads2])
        
        self.logger.info(f"Aligning reads to {reference_fasta}")
        self.logger.debug(f"BWA command: {' '.join(cmd)}")
        
        with open(output_sam, 'w') as sam_file:
            result = subprocess.run(cmd, stdout=sam_file, stderr=subprocess.PIPE, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"BWA alignment failed: {result.stderr}")
        
        if not os.path.isfile(output_sam) or os.path.getsize(output_sam) == 0:
            raise RuntimeError(f"SAM file not created or is empty: {output_sam}")
        
        self.logger.info(f"BWA alignment completed: {output_sam}")
        return output_sam


class SAMtoolsProcessor:
    """Handles SAMtools operations for BAM processing"""
    
    def __init__(self, threads: int = 4, logger: Optional[logging.Logger] = None):
        self.threads = threads
        self.logger = logger or logging.getLogger(__name__)
        self.tool_wrapper = ExternalToolWrapper()
    
    def sam_to_sorted_bam(self, sam_file: str, output_bam: str, 
                         remove_sam: bool = True) -> str:
        """
        Convert SAM to sorted BAM file
        
        Args:
            sam_file: Input SAM file path
            output_bam: Output BAM file path
            remove_sam: Whether to remove SAM file after conversion
            
        Returns:
            Path to sorted BAM file
        """
        if not self.tool_wrapper.validate_tool('samtools'):
            raise RuntimeError("samtools not found in PATH")
        
        temp_bam = f"{output_bam}.tmp"
        
        try:
            self.logger.info(f"Converting {sam_file} to BAM format")
            cmd_view = [
                'samtools', 'view', 
                '-b',
                '-F', '2048',
                '-G', '12',
                '--threads', str(self.threads),
                sam_file
            ]
            
            with open(temp_bam, 'wb') as bam_file:
                result = subprocess.run(cmd_view, stdout=bam_file, stderr=subprocess.PIPE)
            
            if result.returncode != 0:
                raise RuntimeError(f"SAM to BAM conversion failed: {result.stderr.decode()}")
            
            self.logger.info(f"Sorting BAM file to {output_bam}")
            cmd_sort = [
                'samtools', 'sort',
                '-o', output_bam,
                '--threads', str(self.threads),
                temp_bam
            ]
            
            returncode, stdout, stderr = self.tool_wrapper.run_command(cmd_sort)
            
            if os.path.isfile(temp_bam):
                os.remove(temp_bam)
            
            if remove_sam and os.path.isfile(sam_file):
                os.remove(sam_file)
                self.logger.debug(f"Removed SAM file: {sam_file}")
            
        except Exception as e:
            for temp_file in [temp_bam]:
                if os.path.isfile(temp_file):
                    os.remove(temp_file)
            raise e
        
        if not os.path.isfile(output_bam) or os.path.getsize(output_bam) == 0:
            raise RuntimeError(f"BAM file not created or is empty: {output_bam}")
        
        self.logger.info(f"Sorted BAM file created: {output_bam}")
        return output_bam
    
    def index_bam(self, bam_file: str) -> str:
        """
        Create index for BAM file
        
        Args:
            bam_file: Path to BAM file
            
        Returns:
            Path to index file (.bai)
        """
        if not os.path.isfile(bam_file):
            raise FileNotFoundError(f"BAM file not found: {bam_file}")
        
        index_file = f"{bam_file}.bai"
        
        if (os.path.isfile(index_file) and 
            os.path.getmtime(index_file) > os.path.getmtime(bam_file)):
            self.logger.debug(f"BAM index already exists: {index_file}")
            return index_file
        
        self.logger.info(f"Creating BAM index for {bam_file}")
        
        cmd = ['samtools', 'index', bam_file]
        returncode, stdout, stderr = self.tool_wrapper.run_command(cmd)
        
        if not os.path.isfile(index_file):
            raise RuntimeError(f"BAM index creation failed: {bam_file}")
        
        self.logger.info(f"BAM index created: {index_file}")
        return index_file


class TEMapper:
    """High-level mapper for TE analysis with consistent naming"""
    
    def __init__(self, threads: int = 4, logger: Optional[logging.Logger] = None):
        self.threads = threads
        self.logger = logger or logging.getLogger(__name__)
        self.bwa_mapper = BWAMapper(threads, logger)
        self.samtools = SAMtoolsProcessor(threads, logger)
    
    def map_reads_to_ir(self, sample_id: str, reads1: str, reads2: str,
                       ir_fasta: str, output_dir: str, 
                       ir_type: str = "IRL") -> str:
        """
        Map reads to IR sequence and return sorted, indexed BAM
        
        This ensures consistent contig naming in the BAM file to match
        the FASTA sequence ID exactly.
        
        Args:
            sample_id: Sample identifier
            reads1: Forward reads FASTQ path
            reads2: Reverse reads FASTQ path
            ir_fasta: IR reference FASTA path
            output_dir: Output directory for BAM files
            ir_type: Type of IR ('IRL' or 'IRR')
            
        Returns:
            Path to sorted, indexed BAM file
        """
        os.makedirs(output_dir, exist_ok=True)
        
        sam_file = os.path.join(output_dir, f"{sample_id}_{ir_type}.sam")
        bam_file = os.path.join(output_dir, f"{sample_id}_{ir_type}_sorted.bam")
        
        try:
            self.bwa_mapper.align_reads(
                reference_fasta=ir_fasta,
                reads1=reads1,
                reads2=reads2,
                output_sam=sam_file,
                sample_id=sample_id
            )
            
            self.samtools.sam_to_sorted_bam(
                sam_file=sam_file,
                output_bam=bam_file,
                remove_sam=True
            )
            
            self.samtools.index_bam(bam_file)
            
            self.logger.info(f"Successfully mapped {sample_id} to {ir_type}: {bam_file}")
            return bam_file
            
        except Exception as e:
            for temp_file in [sam_file, bam_file, f"{bam_file}.bai"]:
                if os.path.isfile(temp_file):
                    os.remove(temp_file)
                    self.logger.debug(f"Cleaned up failed file: {temp_file}")
            raise e
    
    def map_sample_to_te(self, sample_id: str, reads1: str, reads2: str,
                        irl_fasta: str, irr_fasta: str, 
                        output_dir: str) -> Tuple[str, str]:
        """
        Map a sample to both IRL and IRR references
        
        Args:
            sample_id: Sample identifier
            reads1: Forward reads FASTQ path
            reads2: Reverse reads FASTQ path
            irl_fasta: IRL reference FASTA path
            irr_fasta: IRR reference FASTA path
            output_dir: Output directory for BAM files
            
        Returns:
            Tuple of (irl_bam_path, irr_bam_path)
        """
        self.logger.info(f"Mapping sample {sample_id} to TE references")
        
        irl_bam = self.map_reads_to_ir(
            sample_id=sample_id,
            reads1=reads1,
            reads2=reads2,
            ir_fasta=irl_fasta,
            output_dir=output_dir,
            ir_type="IRL"
        )
        
        irr_bam = self.map_reads_to_ir(
            sample_id=sample_id,
            reads1=reads1,
            reads2=reads2,
            ir_fasta=irr_fasta,
            output_dir=output_dir,
            ir_type="IRR"
        )
        
        return irl_bam, irr_bam


class ParallelMapper:
    """Handles parallel mapping of multiple samples"""
    
    def __init__(self, max_workers: int = 4, threads_per_worker: int = 4,
                 logger: Optional[logging.Logger] = None):
        self.max_workers = max_workers
        self.threads_per_worker = threads_per_worker
        self.logger = logger or logging.getLogger(__name__)
    
    def map_samples_parallel(self, samples: List[Dict], irl_fasta: str, 
                           irr_fasta: str, output_base_dir: str) -> Dict[str, Tuple[str, str]]:
        """
        Map multiple samples in parallel
        
        Args:
            samples: List of sample dictionaries with 'sample_id', 'reads1', 'reads2'
            irl_fasta: IRL reference FASTA path
            irr_fasta: IRR reference FASTA path
            output_base_dir: Base output directory
            
        Returns:
            Dictionary mapping sample_id to (irl_bam, irr_bam) paths
        """
        self.logger.info(f"Starting parallel mapping of {len(samples)} samples with {self.max_workers} workers")
        
        results = {}
        failed_samples = []
        
        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_sample = {}
            for sample in samples:
                sample_id = sample['sample_id']
                sample_output_dir = os.path.join(output_base_dir, sample_id)
                
                future = executor.submit(
                    self._map_single_sample,
                    sample_id,
                    sample['reads1'],
                    sample['reads2'],
                    irl_fasta,
                    irr_fasta,
                    sample_output_dir
                )
                future_to_sample[future] = sample_id
            
            for future in as_completed(future_to_sample):
                sample_id = future_to_sample[future]
                try:
                    irl_bam, irr_bam = future.result()
                    results[sample_id] = (irl_bam, irr_bam)
                    self.logger.info(f"Completed mapping for {sample_id}")
                except Exception as e:
                    self.logger.error(f"Failed to map sample {sample_id}: {str(e)}")
                    failed_samples.append(sample_id)
        
        if failed_samples:
            self.logger.warning(f"Failed to map {len(failed_samples)} samples: {failed_samples}")
        
        self.logger.info(f"Parallel mapping completed. Success: {len(results)}, Failed: {len(failed_samples)}")
        return results
    
    def _map_single_sample(self, sample_id: str, reads1: str, reads2: str,
                          irl_fasta: str, irr_fasta: str, output_dir: str) -> Tuple[str, str]:
        """Map a single sample (used by parallel executor)"""
        mapper = TEMapper(threads=self.threads_per_worker, logger=self.logger)
        return mapper.map_sample_to_te(
            sample_id=sample_id,
            reads1=reads1,
            reads2=reads2,
            irl_fasta=irl_fasta,
            irr_fasta=irr_fasta,
            output_dir=output_dir
        )



if __name__ == "__main__":
    import logging
    from pathlib import Path
    import sys
    import os
    from datetime import datetime
    
    script_dir = Path(__file__).resolve().parent
    project_root = script_dir.parent
    test_output_dir = project_root / "test" / "mapping_test_output"
    
    test_output_dir.mkdir(parents=True, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    test_run_dir = test_output_dir / f"run_{timestamp}"
    test_run_dir.mkdir(exist_ok=True)
    
    log_file = test_run_dir / "mapping_test.log"
    
    file_handler = logging.FileHandler(log_file)
    console_handler = logging.StreamHandler(sys.stdout)
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)
    
    logging.basicConfig(level=logging.DEBUG, handlers=[file_handler, console_handler])
    logger = logging.getLogger(__name__)
    
    logger.info("="*80)
    logger.info("TETyper 2.0 Mapping Module - Comprehensive Testing")
    logger.info("="*80)
    logger.info(f"Project root: {project_root}")
    logger.info(f"Test output directory: {test_run_dir}")
    logger.info(f"Log file: {log_file}")
    logger.info("="*80)
    
    logger.info("\n[TEST 0] Setting up test directory structure...")
    
    ref_dir = test_run_dir / "references"
    reads_dir = test_run_dir / "reads"
    mapping_dir = test_run_dir / "mappings"
    bam_dir = test_run_dir / "bam_files"
    
    for directory in [ref_dir, reads_dir, mapping_dir, bam_dir]:
        directory.mkdir(exist_ok=True)
        logger.info(f"✓ Created directory: {directory.relative_to(project_root)}")
    
    logger.info("\n[TEST 1] Checking external tool availability...")
    
    wrapper = ExternalToolWrapper()
    required_tools = ['bwa', 'samtools']
    tools_available = {}
    tool_versions = {}
    
    for tool in required_tools:
        available = wrapper.validate_tool(tool)
        tools_available[tool] = available
        
        if available:
            logger.info(f"✓ {tool} is available in PATH")
            
            try:
                if tool == 'bwa':
                    returncode, stdout, stderr = wrapper.run_command(['bwa'], check_return=False)
                    version_line = stderr.split('\n')[2] if stderr else "Unknown"
                    tool_versions[tool] = version_line.strip()
                    logger.info(f"  Version: {version_line.strip()}")
                elif tool == 'samtools':
                    returncode, stdout, stderr = wrapper.run_command(['samtools', '--version'], check_return=False)
                    version_line = stdout.split('\n')[0] if stdout else "Unknown"
                    tool_versions[tool] = version_line.strip()
                    logger.info(f"  Version: {version_line.strip()}")
            except Exception as e:
                logger.warning(f"  Could not get version: {e}")
                tool_versions[tool] = "Unknown"
        else:
            logger.error(f"✗ {tool} is NOT available in PATH")
            tool_versions[tool] = "Not Available"
    
    tool_report_file = test_run_dir / "tool_availability.txt"
    with open(tool_report_file, 'w') as f:
        f.write("Tool Availability Report\n")
        f.write("=" * 50 + "\n\n")
        for tool in required_tools:
            status = "Available" if tools_available[tool] else "Not Available"
            f.write(f"{tool}: {status}\n")
            f.write(f"  Version: {tool_versions[tool]}\n\n")
    
    logger.info(f"✓ Tool availability report saved to: {tool_report_file.relative_to(project_root)}")
    
    if not all(tools_available.values()):
        logger.error("Not all required tools are available. Some tests will be skipped.")
    
    logger.info("\n[TEST 2] Creating test reference sequences...")
    
    irl_fasta = ref_dir / "Tn2_IRL_test.fasta"
    irl_sequence = "GGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGG"
    
    with open(irl_fasta, 'w') as f:
        f.write(">Tn2_IRL\n")
        f.write(f"{irl_sequence}\n")
    
    logger.info(f"✓ Created test IRL FASTA: {irl_fasta.relative_to(project_root)}")
    logger.info(f"  Sequence length: {len(irl_sequence)} bp")
    
    irr_fasta = ref_dir / "Tn2_IRR_test.fasta"
    irr_sequence = "GGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGCA"
    
    with open(irr_fasta, 'w') as f:
        f.write(">Tn2_IRR\n")
        f.write(f"{irr_sequence}\n")
    
    logger.info(f"✓ Created test IRR FASTA: {irr_fasta.relative_to(project_root)}")
    logger.info(f"  Sequence length: {len(irr_sequence)} bp")
    
    if tools_available.get('bwa', False):
        logger.info("\n[TEST 3] Testing BWA index creation...")
        
        try:
            indexer = BWAIndexer(logger)
            
            logger.info(f"Creating BWA index for IRL reference...")
            success_irl = indexer.create_index(str(irl_fasta))
            
            if success_irl:
                logger.info("✓ IRL index creation successful")
                
                index_extensions = ['amb', 'ann', 'bwt', 'pac', 'sa']
                index_report = []
                for ext in index_extensions:
                    index_file = f"{irl_fasta}.{ext}"
                    if os.path.isfile(index_file):
                        size = os.path.getsize(index_file)
                        logger.debug(f"  ✓ {ext} file exists ({size} bytes)")
                        index_report.append(f"{ext}: {size} bytes")
                    else:
                        logger.error(f"  ✗ {ext} file missing")
                        index_report.append(f"{ext}: MISSING")
                
                index_report_file = test_run_dir / "index_files_IRL.txt"
                with open(index_report_file, 'w') as f:
                    f.write("BWA Index Files (IRL)\n")
                    f.write("=" * 50 + "\n\n")
                    f.write("\n".join(index_report))
                
                logger.info(f"✓ Index report saved to: {index_report_file.relative_to(project_root)}")
                
                logger.info("Testing index detection (should skip re-indexing)...")
                success_reindex = indexer.create_index(str(irl_fasta))
                if success_reindex:
                    logger.info("✓ Correctly detected existing index")
            
            logger.info(f"Creating BWA index for IRR reference...")
            success_irr = indexer.create_index(str(irr_fasta))
            
            if success_irr:
                logger.info("✓ IRR index creation successful")
            
        except Exception as e:
            logger.error(f"✗ BWA indexing test failed: {str(e)}")
    else:
        logger.warning("[TEST 3] SKIPPED - BWA not available")
    
    logger.info("\n[TEST 4] Creating mock FASTQ files for alignment testing...")
    
    reads1_file = reads_dir / "test_sample_R1.fastq"
    reads2_file = reads_dir / "test_sample_R2.fastq"
    
    mock_reads_r1 = [
        "@READ1",
        "GGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGAAAACCCCGGGG",
        "+",
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        "@READ2",
        "TTTTGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGGGGGCCCC",
        "+",
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        "@READ3",
        "AAAAGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGTTTTAAAA",
        "+",
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
    ]
    
    mock_reads_r2 = [
        "@READ1",
        "CCCCGGGGTTTTCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCC",
        "+",
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        "@READ2",
        "GGGGCCCCCCCCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCAAAA",
        "+",
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        "@READ3",
        "TTTTAAAATTTTCCCTTAACGTGAGTTTTCGTTCCACTGAGCGTCAGACCCCTTTT",
        "+",
        "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
    ]
    
    with open(reads1_file, 'w') as f:
        f.write('\n'.join(mock_reads_r1))
    
    with open(reads2_file, 'w') as f:
        f.write('\n'.join(mock_reads_r2))
    
    logger.info(f"✓ Created mock R1 FASTQ: {reads1_file.relative_to(project_root)}")
    logger.info(f"  Number of reads: {len(mock_reads_r1) // 4}")
    logger.info(f"✓ Created mock R2 FASTQ: {reads2_file.relative_to(project_root)}")
    logger.info(f"  Number of reads: {len(mock_reads_r2) // 4}")
    
    if tools_available.get('bwa', False):
        logger.info("\n[TEST 5] Testing BWA alignment...")
        
        try:
            mapper = BWAMapper(threads=2, logger=logger)
            
            output_sam = mapping_dir / "test_sample_IRL.sam"
            
            logger.info("Running BWA MEM alignment to IRL reference...")
            result_sam = mapper.align_reads(
                reference_fasta=str(irl_fasta),
                reads1=str(reads1_file),
                reads2=str(reads2_file),
                output_sam=str(output_sam),
                sample_id="TEST001"
            )
            
            if os.path.isfile(result_sam):
                sam_size = os.path.getsize(result_sam)
                logger.info(f"✓ SAM file created: {Path(result_sam).relative_to(project_root)}")
                logger.info(f"  File size: {sam_size} bytes")
                
                with open(result_sam, 'r') as f:
                    sam_lines = f.readlines()
                    header_lines = sum(1 for line in sam_lines if line.startswith('@'))
                    alignment_lines = len(sam_lines) - header_lines
                    logger.info(f"  Header lines: {header_lines}")
                    logger.info(f"  Alignment records: {alignment_lines}")
                    
                    alignment_summary_file = test_run_dir / "alignment_summary.txt"
                    with open(alignment_summary_file, 'w') as summary:
                        summary.write("BWA Alignment Summary\n")
                        summary.write("=" * 50 + "\n\n")
                        summary.write(f"SAM file: {result_sam}\n")
                        summary.write(f"File size: {sam_size} bytes\n")
                        summary.write(f"Header lines: {header_lines}\n")
                        summary.write(f"Alignment records: {alignment_lines}\n\n")
                        
                        summary.write("Sample alignments:\n")
                        summary.write("-" * 50 + "\n")
                        count = 0
                        for line in sam_lines:
                            if not line.startswith('@'):
                                summary.write(line)
                                count += 1
                                if count >= 5:
                                    break
                    
                    logger.info(f"✓ Alignment summary saved to: {alignment_summary_file.relative_to(project_root)}")
            else:
                logger.error("✗ SAM file was not created")
            
        except Exception as e:
            logger.error(f"✗ BWA alignment test failed: {str(e)}")
    else:
        logger.warning("[TEST 5] SKIPPED - BWA not available")
    
    if tools_available.get('samtools', False) and tools_available.get('bwa', False):
        logger.info("\n[TEST 6] Testing SAMtools BAM processing...")
        
        try:
            samtools_processor = SAMtoolsProcessor(threads=2, logger=logger)
            
            output_bam = bam_dir / "test_sample_IRL_sorted.bam"
            
            logger.info("Converting SAM to sorted BAM...")
            result_bam = samtools_processor.sam_to_sorted_bam(
                sam_file=str(output_sam),
                output_bam=str(output_bam),
                remove_sam=False
            )
            
            if os.path.isfile(result_bam):
                bam_size = os.path.getsize(result_bam)
                logger.info(f"✓ BAM file created: {Path(result_bam).relative_to(project_root)}")
                logger.info(f"  File size: {bam_size} bytes")
                
                logger.info("Creating BAM index...")
                index_file = samtools_processor.index_bam(str(output_bam))
                
                if os.path.isfile(index_file):
                    index_size = os.path.getsize(index_file)
                    logger.info(f"✓ BAM index created: {Path(index_file).relative_to(project_root)}")
                    logger.info(f"  Index size: {index_size} bytes")
                    
                    logger.info("Testing index detection (should skip re-indexing)...")
                    index_file2 = samtools_processor.index_bam(str(output_bam))
                    if index_file2 == index_file:
                        logger.info("✓ Correctly detected existing index")
                else:
                    logger.error("✗ BAM index file was not created")
                
                logger.info("Retrieving alignment statistics with samtools flagstat...")
                try:
                    cmd = ['samtools', 'flagstat', str(output_bam)]
                    returncode, stdout, stderr = wrapper.run_command(cmd)
                    
                    flagstat_file = test_run_dir / "flagstat_output.txt"
                    with open(flagstat_file, 'w') as f:
                        f.write("Samtools Flagstat Output\n")
                        f.write("=" * 50 + "\n\n")
                        f.write(stdout)
                    
                    logger.info(f"✓ Flagstat output saved to: {flagstat_file.relative_to(project_root)}")
                    logger.info("Alignment statistics:")
                    for line in stdout.strip().split('\n'):
                        logger.info(f"  {line}")
                except Exception as e:
                    logger.warning(f"Could not get flagstat: {e}")
            else:
                logger.error("✗ BAM file was not created")
            
        except Exception as e:
            logger.error(f"✗ SAMtools processing test failed: {str(e)}")
    else:
        logger.warning("[TEST 6] SKIPPED - samtools or BWA not available")
    
    if all(tools_available.values()):
        logger.info("\n[TEST 7] Testing TEMapper high-level interface...")
        
        try:
            te_mapper = TEMapper(threads=2, logger=logger)
            
            logger.info("Testing map_reads_to_ir for IRL...")
            output_dir_mapper = bam_dir / "temapper_output"
            
            irl_bam_result = te_mapper.map_reads_to_ir(
                sample_id="TEST001",
                reads1=str(reads1_file),
                reads2=str(reads2_file),
                ir_fasta=str(irl_fasta),
                output_dir=str(output_dir_mapper),
                ir_type="IRL"
            )
            
            if os.path.isfile(irl_bam_result):
                logger.info(f"✓ IRL BAM created via TEMapper: {Path(irl_bam_result).relative_to(project_root)}")
                
                index_path = f"{irl_bam_result}.bai"
                if os.path.isfile(index_path):
                    logger.info(f"✓ BAM index exists: {Path(index_path).relative_to(project_root)}")
                else:
                    logger.error(f"✗ BAM index missing: {index_path}")
            else:
                logger.error("✗ IRL BAM was not created by TEMapper")
            
            logger.info("Testing map_sample_to_te (IRL + IRR)...")
            output_dir_full = bam_dir / "full_mapping_output"
            
            irl_bam, irr_bam = te_mapper.map_sample_to_te(
                sample_id="TEST001",
                reads1=str(reads1_file),
                reads2=str(reads2_file),
                irl_fasta=str(irl_fasta),
                irr_fasta=str(irr_fasta),
                output_dir=str(output_dir_full)
            )
            
            logger.info(f"✓ Full TE mapping completed")
            logger.info(f"  IRL BAM: {Path(irl_bam).relative_to(project_root)} (exists: {os.path.isfile(irl_bam)})")
            logger.info(f"  IRR BAM: {Path(irr_bam).relative_to(project_root)} (exists: {os.path.isfile(irr_bam)})")
            
            mapping_summary_file = test_run_dir / "temapper_summary.txt"
            with open(mapping_summary_file, 'w') as f:
                f.write("TEMapper Mapping Summary\n")
                f.write("=" * 50 + "\n\n")
                f.write(f"Sample ID: TEST001\n")
                f.write(f"IRL BAM: {irl_bam}\n")
                f.write(f"  Exists: {os.path.isfile(irl_bam)}\n")
                f.write(f"  Size: {os.path.getsize(irl_bam) if os.path.isfile(irl_bam) else 0} bytes\n")
                f.write(f"IRR BAM: {irr_bam}\n")
                f.write(f"  Exists: {os.path.isfile(irr_bam)}\n")
                f.write(f"  Size: {os.path.getsize(irr_bam) if os.path.isfile(irr_bam) else 0} bytes\n")
            
            logger.info(f"✓ Mapping summary saved to: {mapping_summary_file.relative_to(project_root)}")
            
        except Exception as e:
            logger.error(f"✗ TEMapper test failed: {str(e)}")
    else:
        logger.warning("[TEST 7] SKIPPED - Not all tools available")
    
    if all(tools_available.values()):
        logger.info("\n[TEST 8] Testing ParallelMapper with multiple samples...")
        
        try:
            samples = [
                {
                    'sample_id': 'TEST001',
                    'reads1': str(reads1_file),
                    'reads2': str(reads2_file)
                },
                {
                    'sample_id': 'TEST002',
                    'reads1': str(reads1_file),
                    'reads2': str(reads2_file)
                },
            ]
            
            parallel_mapper = ParallelMapper(
                max_workers=2,
                threads_per_worker=1,
                logger=logger
            )
            
            output_dir_parallel = bam_dir / "parallel_output"
            
            logger.info(f"Mapping {len(samples)} samples in parallel...")
            results = parallel_mapper.map_samples_parallel(
                samples=samples,
                irl_fasta=str(irl_fasta),
                irr_fasta=str(irr_fasta),
                output_base_dir=str(output_dir_parallel)
            )
            
            logger.info(f"✓ Parallel mapping completed")
            logger.info(f"  Successfully mapped: {len(results)} samples")
            
            parallel_summary_file = test_run_dir / "parallel_mapping_summary.txt"
            with open(parallel_summary_file, 'w') as f:
                f.write("Parallel Mapping Summary\n")
                f.write("=" * 50 + "\n\n")
                f.write(f"Total samples: {len(samples)}\n")
                f.write(f"Successfully mapped: {len(results)}\n\n")
                
                for sample_id, (irl_bam, irr_bam) in results.items():
                    f.write(f"Sample: {sample_id}\n")
                    f.write(f"  IRL BAM: {irl_bam}\n")
                    f.write(f"    Exists: {os.path.isfile(irl_bam)}\n")
                    f.write(f"  IRR BAM: {irr_bam}\n")
                    f.write(f"    Exists: {os.path.isfile(irr_bam)}\n\n")
                    
                    logger.info(f"  {sample_id}:")
                    logger.info(f"    IRL BAM exists: {os.path.isfile(irl_bam)}")
                    logger.info(f"    IRR BAM exists: {os.path.isfile(irr_bam)}")
            
            logger.info(f"✓ Parallel mapping summary saved to: {parallel_summary_file.relative_to(project_root)}")
            
        except Exception as e:
            logger.error(f"✗ Parallel mapping test failed: {str(e)}")
    else:
        logger.warning("[TEST 8] SKIPPED - Not all tools available")
    
    logger.info("\n" + "="*80)
    logger.info("GENERATING FINAL TEST SUMMARY")
    logger.info("="*80)
    
    summary_file = test_run_dir / "TEST_SUMMARY.txt"
    
    with open(summary_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("TETyper 2.0 Mapping Module - Test Summary\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Test Run: {timestamp}\n")
        f.write(f"Test Directory: {test_run_dir}\n\n")
        
        f.write("Tool Availability:\n")
        f.write("-" * 40 + "\n")
        for tool, available in tools_available.items():
            status = "✓ Available" if available else "✗ Not Available"
            f.write(f"  {tool}: {status}\n")
            f.write(f"    Version: {tool_versions[tool]}\n")
        f.write("\n")
        
        f.write("Test Files Generated:\n")
        f.write("-" * 40 + "\n")
        f.write(f"  References: {ref_dir.relative_to(project_root)}\n")
        f.write(f"  Mock Reads: {reads_dir.relative_to(project_root)}\n")
        f.write(f"  Mappings: {mapping_dir.relative_to(project_root)}\n")
        f.write(f"  BAM Files: {bam_dir.relative_to(project_root)}\n\n")
        
        f.write("Output Files:\n")
        f.write("-" * 40 + "\n")
        for output_file in sorted(test_run_dir.glob("*.txt")):
            if output_file != summary_file:
                f.write(f"  {output_file.name}\n")
        f.write(f"  {log_file.name}\n\n")
        
        f.write("=" * 80 + "\n")
        f.write("All test results have been saved to:\n")
        f.write(f"{test_run_dir}\n")
        f.write("=" * 80 + "\n")
    
    logger.info(f"\n✓ Final test summary saved to: {summary_file.relative_to(project_root)}")
    logger.info("\n" + "="*80)
    logger.info("TEST SUITE COMPLETED SUCCESSFULLY")
    logger.info("="*80)
    logger.info(f"\nAll test outputs saved to: {test_run_dir.relative_to(project_root)}")
    logger.info(f"View complete log: {log_file.relative_to(project_root)}")
    logger.info("="*80)
