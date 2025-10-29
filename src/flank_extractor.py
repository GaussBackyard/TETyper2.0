#!/usr/bin/env python
"""
TETyper 2.0: FlankExtractor Module - Fixed Contig Name Compatibility
Handles mismatches between BAM contig names and FASTA sequence IDs
"""

import pysam
import subprocess
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
from Bio import SeqIO
from Bio.Seq import Seq
import logging
from pathlib import Path


@dataclass 
class ProcessingParameters:
    """Parameters for flanking sequence extraction and quality control"""
    flank_length: int = 10
    min_mapped_length: int = 5
    min_quality: int = 20
    min_reads: int = 10
    min_each_strand: int = 3
    threads: int = 4


class FlankExtractor:
    """
    Enhanced flanking sequence extraction using IR-based detection
    With robust contig name matching between BAM files and FASTA references
    """
    
    def __init__(self, params: ProcessingParameters, logger: Optional[logging.Logger] = None):
        self.params = params
        self.logger = logger or logging.getLogger(__name__)
        
    def _get_bam_contigs(self, bamfile: str) -> List[str]:
        """Get list of all contig names in BAM file"""
        try:
            with pysam.AlignmentFile(bamfile, "rb") as bam:
                return list(bam.references)
        except Exception as e:
            raise ValueError(f"Failed to read BAM file {bamfile}: {str(e)}")
    
    def _find_matching_contig(self, bamfile: str, expected_contig: str) -> str:
        """
        Find the best matching contig name in BAM file for the expected contig
        Handles common naming variations between FASTA ID and BAM contig names
        """
        available_contigs = self._get_bam_contigs(bamfile)
        
        self.logger.debug(f"Available contigs in BAM: {available_contigs}")
        self.logger.debug(f"Looking for contig matching: {expected_contig}")
        
        if expected_contig in available_contigs:
            return expected_contig
            
        variations_to_try = [
            expected_contig.replace("_", ""),
            expected_contig.replace("-", ""),
            expected_contig + "Left",
            expected_contig + "Right",
            expected_contig.replace("IRL", "IRLeft"),
            expected_contig.replace("IRR", "IRRight"),
            expected_contig.upper(),
            expected_contig.lower(),
        ]
        
        for variation in variations_to_try:
            if variation in available_contigs:
                self.logger.info(f"Found matching contig: '{variation}' for expected '{expected_contig}'")
                return variation
                
        for available in available_contigs:
            expected_clean = expected_contig.replace("_", "").replace("-", "").upper()
            available_clean = available.replace("_", "").replace("-", "").upper()
            
            if expected_clean in available_clean or available_clean in expected_clean:
                self.logger.info(f"Found partial match: '{available}' for expected '{expected_contig}'")
                return available
                
        raise ValueError(
            f"No matching contig found for '{expected_contig}' in BAM file {bamfile}.\n"
            f"Available contigs: {available_contigs}\n"
            f"Please ensure the BAM file was aligned against the correct reference sequence."
        )
    
    def extract_flanks_from_bam(self, bamfile: str, ref_contig_fasta: str, 
                               is_reverse_complement: bool = False) -> Dict:
        """
        Extract flanking sequences from BAM alignment with automatic contig name matching
        """
        self.logger.info(f"Extracting flanks from {bamfile} using reference {ref_contig_fasta}")
        
        ref_record, ref_name, ref_length = self._load_reference_sequence(ref_contig_fasta)
        
        self._ensure_bam_index(bamfile)
        
        actual_contig_name = self._find_matching_contig(bamfile, ref_name)
        
        with pysam.AlignmentFile(bamfile, "rb") as samfile:
            left_flanks = self._extract_left_flanks(samfile, actual_contig_name, ref_length)
            right_flanks = self._extract_right_flanks(samfile, actual_contig_name, ref_length)
            
        result = self._process_extracted_flanks(
            left_flanks, right_flanks, is_reverse_complement, ref_name, ref_length
        )
        
        result['bam_contig_used'] = actual_contig_name
        result['fasta_sequence_id'] = ref_name
        
        self.logger.info(f"Extracted {len(result['left'])} left and {len(result['right'])} right flanks")
        return result
    
    def _load_reference_sequence(self, ref_contig_fasta: str) -> Tuple[object, str, int]:
        """Load IR reference from FASTA and return record, name, and length"""
        try:
            records = list(SeqIO.parse(ref_contig_fasta, "fasta"))
            if len(records) != 1:
                raise ValueError(f"Expected exactly one record in {ref_contig_fasta}, found {len(records)}")
            
            ref_record = records[0]
            ref_name = ref_record.id
            ref_length = len(ref_record.seq)
            
            self.logger.debug(f"Loaded reference {ref_name} with length {ref_length} from {ref_contig_fasta}")
            return ref_record, ref_name, ref_length
            
        except Exception as e:
            raise ValueError(f"Failed to load reference from {ref_contig_fasta}: {str(e)}")
    
    def _ensure_bam_index(self, bamfile: str) -> None:
        """Ensure BAM file is indexed, create index if needed"""
        try:
            with pysam.AlignmentFile(bamfile, "rb") as _tmp:
                pass
        except FileNotFoundError as e:
            raise FileNotFoundError(f"BAM file not found: {bamfile}") from e
        
        try:
            pysam.idxstats(bamfile)
        except Exception:
            self.logger.info(f"Creating BAM index for {bamfile}")
            subprocess.run(["samtools", "index", bamfile], check=True)
    
    def _extract_left_flanks(self, samfile: pysam.AlignmentFile, ref_name: str, 
                           ref_length: int) -> List[Tuple[str, int]]:
        """Extract left flanking sequences from reads overlapping the start boundary"""
        left_flanks = []
        ref_start = 1
        
        try:
            for read in samfile.fetch(ref_name, ref_start - 1, ref_start):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                    
                ref_positions = read.get_reference_positions(full_length=True)
                
                if not ref_positions:
                    continue
                    
                if (ref_start - 1) in ref_positions:
                    mapped_downstream = sum(
                        1 for pos in ref_positions 
                        if (pos is not None and pos >= (ref_start - 1 + (self.params.min_mapped_length - 1)))
                    )
                    
                    if mapped_downstream > 0:
                        try:
                            idx = ref_positions.index(ref_start - 1)
                            
                            flank_start = max(0, idx - self.params.flank_length)
                            flank_seq = read.query_sequence[flank_start:idx]
                            flank_qual = read.query_qualities[flank_start:idx] if read.query_qualities else []
                            
                            if (len(flank_seq) == self.params.flank_length and 
                                (not flank_qual or min(flank_qual) >= self.params.min_quality)):
                                
                                strand = -1 if read.is_reverse else 1
                                left_flanks.append((flank_seq, strand))
                                
                        except (ValueError, IndexError):
                            continue
        except Exception as e:
            self.logger.error(f"Error extracting left flanks from contig {ref_name}: {str(e)}")
            raise
        
        return left_flanks
    
    def _extract_right_flanks(self, samfile: pysam.AlignmentFile, ref_name: str, 
                            ref_length: int) -> List[Tuple[str, int]]:
        """Extract right flanking sequences from reads overlapping the end boundary"""
        right_flanks = []
        ref_end = ref_length
        
        try:
            for read in samfile.fetch(ref_name, ref_end - 1, ref_end):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                    
                ref_positions = read.get_reference_positions(full_length=True)
                
                if not ref_positions:
                    continue
                    
                if (ref_end - 1) in ref_positions:
                    mapped_upstream = sum(
                        1 for pos in ref_positions 
                        if (pos is not None and pos <= (ref_end - 1 - (self.params.min_mapped_length - 1)))
                    )
                    
                    if mapped_upstream > 0:
                        try:
                            idx = ref_positions.index(ref_end - 1)
                            
                            flank_end = min(len(read.query_sequence), idx + 1 + self.params.flank_length)
                            flank_seq = read.query_sequence[idx + 1:flank_end]
                            flank_qual = read.query_qualities[idx + 1:flank_end] if read.query_qualities else []
                            
                            if (len(flank_seq) == self.params.flank_length and 
                                (not flank_qual or min(flank_qual) >= self.params.min_quality)):
                                
                                strand = -1 if read.is_reverse else 1
                                right_flanks.append((flank_seq, strand))
                                
                        except (ValueError, IndexError):
                            continue
        except Exception as e:
            self.logger.error(f"Error extracting right flanks from contig {ref_name}: {str(e)}")
            raise
        
        return right_flanks
    
    def _process_extracted_flanks(self, left_flanks: List[Tuple[str, int]], 
                                right_flanks: List[Tuple[str, int]], 
                                is_reverse_complement: bool,
                                ref_name: str, ref_length: int) -> Dict:
        """Process extracted flanks and apply quality control filters"""
        left_counts = Counter([seq for seq, strand in left_flanks])
        right_counts = Counter([seq for seq, strand in right_flanks])
        
        left_strand = defaultdict(lambda: {'forward': 0, 'reverse': 0})
        for seq, strand in left_flanks:
            if strand == 1:
                left_strand[seq]['forward'] += 1
            else:
                left_strand[seq]['reverse'] += 1
        
        right_strand = defaultdict(lambda: {'forward': 0, 'reverse': 0})
        for seq, strand in right_flanks:
            if strand == 1:
                right_strand[seq]['forward'] += 1
            else:
                right_strand[seq]['reverse'] += 1
        
        result = {
            'reference_name': ref_name,
            'reference_length': ref_length,
            'left': [], 
            'right': []
        }
        
        for seq in left_counts:
            count = left_counts[seq]
            forward_reads = left_strand[seq]['forward']
            reverse_reads = left_strand[seq]['reverse']
            
            passes_filter = (count >= self.params.min_reads and 
                           forward_reads >= self.params.min_each_strand and 
                           reverse_reads >= self.params.min_each_strand)
            
            if is_reverse_complement:
                processed_seq = str(Seq(seq).reverse_complement())
                forward_reads, reverse_reads = reverse_reads, forward_reads
            else:
                processed_seq = seq
            
            result['left'].append({
                'flank': processed_seq,
                'count': count,
                'forward_reads': forward_reads,
                'reverse_reads': reverse_reads,
                'passes_filter': passes_filter
            })
        
        for seq in right_counts:
            count = right_counts[seq]
            forward_reads = right_strand[seq]['forward']
            reverse_reads = right_strand[seq]['reverse']
            
            passes_filter = (count >= self.params.min_reads and 
                           forward_reads >= self.params.min_each_strand and 
                           reverse_reads >= self.params.min_each_strand)
            
            if is_reverse_complement:
                processed_seq = str(Seq(seq).reverse_complement())
                forward_reads, reverse_reads = reverse_reads, forward_reads
            else:
                processed_seq = seq
            
            result['right'].append({
                'flank': processed_seq,
                'count': count,
                'forward_reads': forward_reads,
                'reverse_reads': reverse_reads,
                'passes_filter': passes_filter
            })
        
        return result
    
    def extract_flanks_for_te(self, sample_id: str, bam_left: str, bam_right: str,
                             irl_fasta: str, irr_fasta: str, te_name: str) -> Dict:
        """Extract flanks for a specific TE using both IRL and IRR BAM files"""
        results = {
            'sample_id': sample_id,
            'te_name': te_name,
            'left_flanks': [],
            'right_flanks': []
        }
        
        if bam_left and Path(bam_left).exists():
            left_result = self.extract_flanks_from_bam(
                bamfile=bam_left,
                ref_contig_fasta=irl_fasta,
                is_reverse_complement=False
            )
            
            for record in left_result['left']:
                if record['passes_filter']:
                    results['left_flanks'].append({
                        'sample': sample_id,
                        'ref': left_result['reference_name'],
                        'biological_side': 'left',
                        'flank': record['flank'],
                        'count': record['count'],
                        'forward_reads': record['forward_reads'],
                        'reverse_reads': record['reverse_reads'],
                        'passes_filter': record['passes_filter']
                    })
        
        if bam_right and Path(bam_right).exists():
            right_result = self.extract_flanks_from_bam(
                bamfile=bam_right,
                ref_contig_fasta=irr_fasta,
                is_reverse_complement=True
            )
            
            for record in right_result['left']:
                if record['passes_filter']:
                    results['right_flanks'].append({
                        'sample': sample_id,
                        'ref': right_result['reference_name'],
                        'biological_side': 'right',
                        'flank': record['flank'],
                        'count': record['count'],
                        'forward_reads': record['forward_reads'],
                        'reverse_reads': record['reverse_reads'],
                        'passes_filter': record['passes_filter']
                    })
        
        self.logger.info(f"Sample {sample_id}: {len(results['left_flanks'])} left flanks, "
                        f"{len(results['right_flanks'])} right flanks")
        print(results)
        return results


def test_flank_extractor_with_existing_bams():
    """Test function using your existing BAM files with contig name matching"""
    
    params = ProcessingParameters(
        flank_length=10,
        min_mapped_length=5,
        min_quality=20,
        min_reads=3,
        min_each_strand=1
    )
    
    extractor = FlankExtractor(params)
    
    
    sample_results = extractor.extract_flanks_for_te(
        sample_id="CAV1016",
        bam_left="/scratch/jl9gx/TETyper/AllCAV_Tn2IR/CAV1016_Tn2IR_Left_sorted.bam",
        bam_right="/scratch/jl9gx/TETyper/AllCAV_Tn2IR/CAV1016_Tn2IR_Right_sorted.bam", 
        irl_fasta="/home/jl9gx/plasmid/TETyper2.0/Tn2_IRL.fasta",
        irr_fasta="/home/jl9gx/plasmid/TETyper2.0/Tn2_IRR.fasta",
        te_name="Tn2"
    )
    
    print(f"Test completed for {sample_results['sample_id']}")
    print(f"Left flanks found: {len(sample_results['left_flanks'])}")
    print(f"Right flanks found: {len(sample_results['right_flanks'])}")
    
    return sample_results


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    test_flank_extractor_with_existing_bams()