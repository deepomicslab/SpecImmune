#!/usr/bin/env python3
"""
HLA Long-Read RNA-seq Simulator with RNA Degradation

模拟不同测序平台的长读长RNA-seq数据，考虑5'端RNA降解问题

Author: Claude
Date: 2024
"""

import argparse
import random
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional
import gzip

# ============================================================================
# Constants and Configuration
# ============================================================================

CLASS1_GENES = ['A', 'B', 'C', 'E', 'F', 'G']
CLASS1_PSEUDOGENES = ['H', 'J', 'K', 'L', 'N', 'P', 'S', 'T', 'U', 'V', 'W']
CLASS2_GENES = ['DRA', 'DQA1', 'DQA2', 'DQB1', 'DPA1', 'DPA2', 'DPB1', 'DPB2', 
                'DMA', 'DMB', 'DOA', 'DOB', 'DRB1', 'DRB3', 'DRB4']
NON_HLA_GENES = ['MICA', 'MICB', 'TAP1', 'TAP2']

HLA_CLASSES = {
    'class1_genes': CLASS1_GENES,
    'class1_pseudogenes': CLASS1_PSEUDOGENES,
    'class2_genes': CLASS2_GENES,
    'non_hla_genes': NON_HLA_GENES
}

# Sequencing platform error profiles
PLATFORM_PROFILES = {
    'pacbio_clr': {
        'name': 'PacBio CLR',
        'error_rate': 0.13,
        'insertion_rate': 0.06,
        'deletion_rate': 0.04,
        'substitution_rate': 0.03,
        'mean_length_ratio': 0.95,
        'length_std_ratio': 0.15
    },
    'pacbio_hifi': {
        'name': 'PacBio HiFi',
        'error_rate': 0.005,
        'insertion_rate': 0.002,
        'deletion_rate': 0.002,
        'substitution_rate': 0.001,
        'mean_length_ratio': 0.98,
        'length_std_ratio': 0.05
    },
    'ont': {
        'name': 'Oxford Nanopore',
        'error_rate': 0.10,
        'insertion_rate': 0.03,
        'deletion_rate': 0.04,
        'substitution_rate': 0.03,
        'mean_length_ratio': 0.92,
        'length_std_ratio': 0.20
    },
    'ont_r10': {
        'name': 'Oxford Nanopore R10',
        'error_rate': 0.05,
        'insertion_rate': 0.015,
        'deletion_rate': 0.02,
        'substitution_rate': 0.015,
        'mean_length_ratio': 0.95,
        'length_std_ratio': 0.12
    }
}

BASES = ['A', 'T', 'G', 'C']

# ============================================================================
# Data Classes
# ============================================================================

@dataclass
class HLAAllele:
    """HLA allele information"""
    hla_id: str
    allele_name: str
    length: int
    sequence: str
    gene: str
    
    @classmethod
    def from_fasta_header(cls, header: str, sequence: str) -> Optional['HLAAllele']:
        """
        Parse HLA allele from FASTA header
        Format: ">HLA:HLA00001 A*01:01:01:01 1098 bp"
        """
        try:
            parts = header.strip().lstrip('>').split()
            if len(parts) < 4:
                return None
            
            hla_id = parts[0]  # HLA:HLA00001
            allele_name = parts[1]  # A*01:01:01:01
            length = int(parts[2])  # 1098
            
            # Extract gene name from allele name
            gene = allele_name.split('*')[0]
            
            return cls(
                hla_id=hla_id,
                allele_name=allele_name,
                length=length,
                sequence=sequence.upper(),
                gene=gene
            )
        except (IndexError, ValueError) as e:
            print(f"Warning: Failed to parse header: {header}, Error: {e}")
            return None

@dataclass
class SimulatedRead:
    """Simulated sequencing read"""
    read_id: str
    sequence: str
    quality: str
    original_start: int
    original_end: int
    degradation_rate: float
    platform: str
    source_allele: str

# ============================================================================
# FASTA Parsing
# ============================================================================

def parse_fasta(fasta_file: str) -> Dict[str, HLAAllele]:
    """Parse FASTA file and return dictionary of HLA alleles"""
    alleles = {}
    
    open_func = gzip.open if fasta_file.endswith('.gz') else open
    mode = 'rt' if fasta_file.endswith('.gz') else 'r'
    
    with open_func(fasta_file, mode) as f:
        current_header = None
        current_seq = []
        
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            if line.startswith('>'):
                # Save previous sequence
                if current_header and current_seq:
                    seq = ''.join(current_seq)
                    allele = HLAAllele.from_fasta_header(current_header, seq)
                    if allele:
                        alleles[allele.allele_name] = allele
                
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_header and current_seq:
            seq = ''.join(current_seq)
            allele = HLAAllele.from_fasta_header(current_header, seq)
            if allele:
                alleles[allele.allele_name] = allele
    
    return alleles

def classify_alleles(alleles: Dict[str, HLAAllele]) -> Dict[str, List[HLAAllele]]:
    """Classify alleles by HLA class"""
    classified = {
        'class1_genes': [],
        'class1_pseudogenes': [],
        'class2_genes': [],
        'non_hla_genes': []
    }
    
    for allele in alleles.values():
        gene = allele.gene
        
        if gene in CLASS1_GENES:
            classified['class1_genes'].append(allele)
        elif gene in CLASS1_PSEUDOGENES:
            classified['class1_pseudogenes'].append(allele)
        elif gene in CLASS2_GENES:
            classified['class2_genes'].append(allele)
        elif gene in NON_HLA_GENES:
            classified['non_hla_genes'].append(allele)
    
    return classified

# ============================================================================
# RNA Degradation Simulation
# ============================================================================

def apply_5prime_degradation(sequence: str, degradation_rate: float) -> Tuple[str, int, int]:
    """
    Simulate 5' end RNA degradation
    
    RNA degradation typically occurs from the 5' end due to:
    - 5' to 3' exonuclease activity
    - Decapping followed by degradation
    
    Args:
        sequence: Original RNA sequence
        degradation_rate: Rate of 5' degradation (0.0 - 1.0)
    
    Returns:
        Tuple of (degraded_sequence, start_position, end_position)
    """
    seq_len = len(sequence)
    
    if seq_len == 0:
        return sequence, 0, 0
    
    # Calculate maximum degradation length based on rate
    # The degradation follows an exponential distribution
    max_degradation = int(seq_len * degradation_rate)
    
    if max_degradation == 0:
        return sequence, 0, seq_len
    
    # Use exponential distribution for degradation length
    # Lambda = 1 / (degradation_rate * seq_len / 2) to center around expected degradation
    lambda_param = 2.0 / max_degradation if max_degradation > 0 else 1.0
    degradation_length = min(
        int(random.expovariate(lambda_param)),
        max_degradation,
        seq_len - 50  # Keep at least 50bp
    )
    
    # Ensure we keep a reasonable length
    degradation_length = max(0, min(degradation_length, seq_len - 100))
    
    start_pos = degradation_length
    end_pos = seq_len
    
    return sequence[start_pos:end_pos], start_pos, end_pos

def apply_3prime_degradation(sequence: str, degradation_rate: float) -> Tuple[str, int]:
    """
    Simulate minor 3' end degradation (less common than 5')
    
    Returns:
        Tuple of (degraded_sequence, new_end_position)
    """
    seq_len = len(sequence)
    
    # 3' degradation is typically less severe
    effective_rate = degradation_rate * 0.3
    max_degradation = int(seq_len * effective_rate)
    
    if max_degradation == 0:
        return sequence, seq_len
    
    lambda_param = 3.0 / max_degradation if max_degradation > 0 else 1.0
    degradation_length = min(
        int(random.expovariate(lambda_param)),
        max_degradation,
        seq_len - 50
    )
    
    degradation_length = max(0, degradation_length)
    new_end = seq_len - degradation_length
    
    return sequence[:new_end], new_end

# ============================================================================
# Sequencing Error Simulation
# ============================================================================

def introduce_sequencing_errors(sequence: str, platform: str) -> str:
    """
    Introduce platform-specific sequencing errors
    
    Args:
        sequence: Input sequence
        platform: Sequencing platform name
    
    Returns:
        Sequence with introduced errors
    """
    profile = PLATFORM_PROFILES.get(platform, PLATFORM_PROFILES['ont'])
    
    result = []
    i = 0
    
    while i < len(sequence):
        base = sequence[i]
        rand = random.random()
        
        # Insertion
        if rand < profile['insertion_rate']:
            result.append(base)
            result.append(random.choice(BASES))
            i += 1
        # Deletion
        elif rand < profile['insertion_rate'] + profile['deletion_rate']:
            i += 1
        # Substitution
        elif rand < profile['error_rate']:
            other_bases = [b for b in BASES if b != base]
            result.append(random.choice(other_bases))
            i += 1
        else:
            result.append(base)
            i += 1
    
    return ''.join(result)

def generate_quality_string(length: int, platform: str) -> str:
    """
    Generate quality string for simulated reads
    
    Args:
        length: Length of the read
        platform: Sequencing platform
    
    Returns:
        Quality string in Phred+33 format
    """
    profile = PLATFORM_PROFILES.get(platform, PLATFORM_PROFILES['ont'])
    
    # Base quality depends on error rate
    # Q = -10 * log10(error_rate)
    base_error = profile['error_rate']
    
    if base_error > 0:
        import math
        mean_q = -10 * math.log10(base_error)
    else:
        mean_q = 40  # High quality for very low error rates
    
    mean_q = min(max(mean_q, 5), 40)  # Clamp to reasonable range
    
    qualities = []
    for _ in range(length):
        # Add some variation
        q = int(random.gauss(mean_q, 3))
        q = min(max(q, 0), 40)
        qualities.append(chr(q + 33))
    
    return ''.join(qualities)

# ============================================================================
# Read Simulation
# ============================================================================

def simulate_reads(
    allele: HLAAllele,
    num_reads: int,
    degradation_rate: float,
    platform: str,
    read_prefix: str
) -> List[SimulatedRead]:
    """
    Simulate long reads from an HLA allele
    
    Args:
        allele: HLA allele to simulate from
        num_reads: Number of reads to generate
        degradation_rate: 5' degradation rate
        platform: Sequencing platform
        read_prefix: Prefix for read IDs
    
    Returns:
        List of simulated reads
    """
    reads = []
    sequence = allele.sequence
    
    for i in range(num_reads):
        # Apply 5' degradation
        degraded_seq, start_pos, end_pos = apply_5prime_degradation(
            sequence, degradation_rate
        )
        
        # Optionally apply 3' degradation
        if random.random() < 0.3:  # 30% chance of 3' degradation
            degraded_seq, relative_end = apply_3prime_degradation(
                degraded_seq, degradation_rate
            )
            end_pos = start_pos + relative_end
        
        if len(degraded_seq) < 50:
            continue
        
        # Introduce sequencing errors
        error_seq = introduce_sequencing_errors(degraded_seq, platform)
        
        # Generate quality string
        quality = generate_quality_string(len(error_seq), platform)
        
        read_id = f"{read_prefix}_{i+1}_deg{degradation_rate:.1f}_{platform}"
        
        read = SimulatedRead(
            read_id=read_id,
            sequence=error_seq,
            quality=quality,
            original_start=start_pos,
            original_end=end_pos,
            degradation_rate=degradation_rate,
            platform=platform,
            source_allele=allele.allele_name
        )
        reads.append(read)
    
    return reads

# ============================================================================
# Output Functions
# ============================================================================

def write_fasta(reads: List[SimulatedRead], output_file: str):
    """Write reads to FASTA file"""
    with open(output_file, 'w') as f:
        for read in reads:
            f.write(f">{read.read_id} source={read.source_allele} "
                   f"start={read.original_start} end={read.original_end}\n")
            # Write sequence in 80-character lines
            seq = read.sequence
            for i in range(0, len(seq), 80):
                f.write(seq[i:i+80] + '\n')

def write_fastq(reads: List[SimulatedRead], output_file: str):
    """Write reads to FASTQ file"""
    with open(output_file, 'w') as f:
        for read in reads:
            f.write(f"@{read.read_id} source={read.source_allele} "
                   f"start={read.original_start} end={read.original_end}\n")
            f.write(f"{read.sequence}\n")
            f.write("+\n")
            f.write(f"{read.quality}\n")

def write_truth_allele_fasta(allele: HLAAllele, output_file: str):
    """Write truth allele sequence to FASTA file"""
    with open(output_file, 'w') as f:
        f.write(f">{allele.hla_id} {allele.allele_name} {allele.length} bp\n")
        seq = allele.sequence
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + '\n')

def write_truth_record(
    selected_alleles: Dict[str, List[HLAAllele]],
    output_file: str
):
    """Write truth record TSV file"""
    with open(output_file, 'w') as f:
        # Write header
        f.write("allele_name\tallele_id\tgene\tallele_class\tsequence_length\n")
        
        for allele_class, alleles in selected_alleles.items():
            for allele in alleles:
                f.write(f"{allele.allele_name}\t{allele.hla_id}\t"
                       f"{allele.gene}\t{allele_class}\t{allele.length}\n")

# ============================================================================
# Main Simulation Pipeline
# ============================================================================

def run_simulation(
    fasta_file: str,
    output_dir: str,
    alleles_per_class: int = 200,
    reads_per_allele: int = 100,
    degradation_rates: List[float] = None,
    platforms: List[str] = None,
    seed: int = None,
    output_format: str = 'fastq'
):
    """
    Run the complete simulation pipeline
    
    Args:
        fasta_file: Input HLA FASTA file
        output_dir: Output directory
        alleles_per_class: Number of alleles to sample per class
        reads_per_allele: Number of reads per allele
        degradation_rates: List of degradation rates to simulate
        platforms: List of platforms to simulate
        seed: Random seed for reproducibility
        output_format: Output format ('fasta' or 'fastq')
    """
    if seed is not None:
        random.seed(seed)
    
    if degradation_rates is None:
        degradation_rates = [0.1, 0.2, 0.3, 0.4, 0.5]
    
    if platforms is None:
        platforms = list(PLATFORM_PROFILES.keys())
    
    # Create output directories
    os.makedirs(output_dir, exist_ok=True)
    truth_dir = os.path.join(output_dir, 'truth_sequences')
    reads_dir = os.path.join(output_dir, 'simulated_reads')
    os.makedirs(truth_dir, exist_ok=True)
    os.makedirs(reads_dir, exist_ok=True)
    
    print(f"Parsing FASTA file: {fasta_file}")
    alleles = parse_fasta(fasta_file)
    print(f"Parsed {len(alleles)} alleles")
    
    # Classify alleles
    classified = classify_alleles(alleles)
    
    print("\nAllele counts by class:")
    for class_name, class_alleles in classified.items():
        print(f"  {class_name}: {len(class_alleles)} alleles")
    
    # Sample alleles from each class
    selected_alleles = {}
    
    for class_name, class_alleles in classified.items():
        if len(class_alleles) == 0:
            print(f"Warning: No alleles found for {class_name}")
            selected_alleles[class_name] = []
            continue
        
        n_sample = min(alleles_per_class, len(class_alleles))
        selected = random.sample(class_alleles, n_sample)
        selected_alleles[class_name] = selected
        print(f"Selected {n_sample} alleles from {class_name}")
    
    # Write truth record
    truth_record_file = os.path.join(output_dir, 'truth_record.tsv')
    write_truth_record(selected_alleles, truth_record_file)
    print(f"\nWrote truth record to: {truth_record_file}")
    
    # Write truth sequences and simulate reads
    total_reads = 0
    
    for class_name, alleles_list in selected_alleles.items():
        class_truth_dir = os.path.join(truth_dir, class_name)
        os.makedirs(class_truth_dir, exist_ok=True)
        
        for allele in alleles_list:
            # Write truth sequence
            safe_name = allele.allele_name.replace('*', '_').replace(':', '_')
            truth_file = os.path.join(class_truth_dir, f"{safe_name}.fasta")
            write_truth_allele_fasta(allele, truth_file)
            
            # Simulate reads for each platform and degradation rate
            for platform in platforms:
                platform_dir = os.path.join(reads_dir, platform)
                os.makedirs(platform_dir, exist_ok=True)
                
                for deg_rate in degradation_rates:
                    deg_dir = os.path.join(platform_dir, f"deg_{deg_rate:.1f}")
                    os.makedirs(deg_dir, exist_ok=True)
                    
                    read_prefix = f"{safe_name}"
                    reads = simulate_reads(
                        allele=allele,
                        num_reads=reads_per_allele,
                        degradation_rate=deg_rate,
                        platform=platform,
                        read_prefix=read_prefix
                    )
                    
                    if reads:
                        ext = 'fq' if output_format == 'fastq' else 'fa'
                        output_file = os.path.join(deg_dir, f"{safe_name}.{ext}")
                        
                        if output_format == 'fastq':
                            write_fastq(reads, output_file)
                        else:
                            write_fasta(reads, output_file)
                        
                        total_reads += len(reads)
    
    print(f"\nSimulation complete!")
    print(f"Total reads generated: {total_reads}")
    print(f"Output directory: {output_dir}")
    
    # Write summary statistics
    summary_file = os.path.join(output_dir, 'simulation_summary.txt')
    with open(summary_file, 'w') as f:
        f.write("HLA Long-Read RNA-seq Simulation Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Input file: {fasta_file}\n")
        f.write(f"Total alleles parsed: {len(alleles)}\n")
        f.write(f"Alleles per class (max): {alleles_per_class}\n")
        f.write(f"Reads per allele: {reads_per_allele}\n")
        f.write(f"Degradation rates: {degradation_rates}\n")
        f.write(f"Platforms: {platforms}\n")
        f.write(f"Random seed: {seed}\n")
        f.write(f"Output format: {output_format}\n\n")
        
        f.write("Selected alleles by class:\n")
        for class_name, alleles_list in selected_alleles.items():
            f.write(f"  {class_name}: {len(alleles_list)}\n")
        
        f.write(f"\nTotal reads generated: {total_reads}\n")
    
    print(f"Summary written to: {summary_file}")

# ============================================================================
# Command Line Interface
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Simulate long-read RNA-seq data from HLA sequences with RNA degradation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python simulate_hla_reads.py -i hla_nuc.fasta -o output_dir

  # With specific platforms and degradation rates
  python simulate_hla_reads.py -i hla_nuc.fasta -o output_dir \\
      --platforms pacbio_hifi ont_r10 \\
      --degradation-rates 0.1 0.2 0.3

  # Customize number of alleles and reads
  python simulate_hla_reads.py -i hla_nuc.fasta -o output_dir \\
      --alleles-per-class 100 \\
      --reads-per-allele 50

HLA Classes:
  class1_genes: A, B, C, E, F, G
  class1_pseudogenes: H, J, K, L, N, P, S, T, U, V, W
  class2_genes: DRA, DQA1, DQA2, DQB1, DPA1, DPA2, DPB1, DPB2, DMA, DMB, DOA, DOB, DRB1, DRB3, DRB4
  non_hla_genes: MICA, MICB, TAP1, TAP2

Platforms:
  pacbio_clr: PacBio CLR (high error rate ~13%)
  pacbio_hifi: PacBio HiFi (low error rate ~0.5%)
  ont: Oxford Nanopore (error rate ~10%)
  ont_r10: Oxford Nanopore R10 (error rate ~5%)
        """
    )
    
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='Input HLA FASTA file (hla_nuc.fasta)'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output directory'
    )
    
    parser.add_argument(
        '--alleles-per-class',
        type=int,
        default=200,
        help='Number of alleles to sample per HLA class (default: 200)'
    )
    
    parser.add_argument(
        '--reads-per-allele',
        type=int,
        default=100,
        help='Number of reads to generate per allele (default: 100)'
    )
    
    parser.add_argument(
        '--degradation-rates',
        type=float,
        nargs='+',
        default=[0.1, 0.2, 0.3, 0.4, 0.5],
        help='5\' degradation rates to simulate (default: 0.1 0.2 0.3 0.4 0.5)'
    )
    
    parser.add_argument(
        '--platforms',
        nargs='+',
        choices=list(PLATFORM_PROFILES.keys()),
        default=list(PLATFORM_PROFILES.keys()),
        help=f'Sequencing platforms to simulate (default: all)'
    )
    
    parser.add_argument(
        '--seed',
        type=int,
        default=42,
        help='Random seed for reproducibility (default: 42)'
    )
    
    parser.add_argument(
        '--format',
        choices=['fasta', 'fastq'],
        default='fastq',
        help='Output format (default: fastq)'
    )
    
    args = parser.parse_args()
    
    # Validate input file
    if not os.path.exists(args.input):
        print(f"Error: Input file not found: {args.input}")
        sys.exit(1)
    
    # Validate degradation rates
    for rate in args.degradation_rates:
        if not 0.0 <= rate <= 1.0:
            print(f"Error: Degradation rate must be between 0.0 and 1.0: {rate}")
            sys.exit(1)
    
    # Run simulation
    run_simulation(
        fasta_file=args.input,
        output_dir=args.output,
        alleles_per_class=args.alleles_per_class,
        reads_per_allele=args.reads_per_allele,
        degradation_rates=args.degradation_rates,
        platforms=args.platforms,
        seed=args.seed,
        output_format=args.format
    )

if __name__ == '__main__':
    main()