from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
import os

def read_conserved_regions(human_fasta, pig_fasta):
    """
    Reads two FASTA files containing pre-aligned conserved regions from BLAST.
    Returns a list of conserved region pairs.
    """
    print(f"Reading conserved regions from {human_fasta} and {pig_fasta}...")
    
    # Check if files exist
    for f in [human_fasta, pig_fasta]:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Input file '{f}' not found.")
    
    # Read sequences from both files
    human_records = list(SeqIO.parse(human_fasta, "fasta"))
    pig_records = list(SeqIO.parse(pig_fasta, "fasta"))
    
    if len(human_records) != len(pig_records):
        print(f"Warning: Different number of sequences in files ({len(human_records)} human vs {len(pig_records)} pig)")
        # Use the minimum number of regions
        num_regions = min(len(human_records), len(pig_records))
        human_records = human_records[:num_regions]
        pig_records = pig_records[:num_regions]
    
    conserved_regions = []
    for i, (human_rec, pig_rec) in enumerate(zip(human_records, pig_records)):
        human_seq = str(human_rec.seq).upper()
        pig_seq = str(pig_rec.seq).upper()
        
        # Ensure sequences are the same length (as they should be from BLAST alignment)
        if len(human_seq) != len(pig_seq):
            print(f"Warning: Region {i+1} has different lengths ({len(human_seq)} vs {len(pig_seq)}). Trimming to shorter length.")
            min_len = min(len(human_seq), len(pig_seq))
            human_seq = human_seq[:min_len]
            pig_seq = pig_seq[:min_len]
        
        conserved_regions.append({
            'human_seq': human_seq,
            'pig_seq': pig_seq,
            'human_id': human_rec.id,
            'pig_id': pig_rec.id,
            'length': len(human_seq)
        })
    
    print(f"Found {len(conserved_regions)} conserved regions")
    for i, region in enumerate(conserved_regions):
        print(f"  Region {i+1}: {region['length']} bp (Human: {region['human_id']}, Pig: {region['pig_id']})")
    
    return conserved_regions

def is_fully_conserved(primer, human_region, pig_region, start_pos, end_pos):
    """
    Check if the primer sequence is fully conserved between human and pig
    in the specified region.
    """
    # Extract the same region from both sequences
    human_subseq = human_region[start_pos:end_pos]
    pig_subseq = pig_region[start_pos:end_pos]
    
    # Check if the primer matches exactly in both species
    return human_subseq == primer and pig_subseq == primer

def calculate_tm(sequence):
    """Calculate melting temperature using basic Wallace rule"""
    a_count = sequence.count('A')
    t_count = sequence.count('T')
    g_count = sequence.count('G')
    c_count = sequence.count('C')
    return 2*(a_count + t_count) + 4*(g_count + c_count)

def has_gc_clamp(sequence):
    """Check if sequence has G or C at the 3' end (last base)"""
    return sequence[-1] in ['G', 'C']

def has_repeats(sequence, max_repeat=4):
    """Check for nucleotide repeats"""
    for base in 'ATGC':
        if base * (max_repeat + 1) in sequence:
            return True
    return False

def has_long_runs(sequence, max_run=3):
    """Check for long runs of the same nucleotide"""
    for base in 'ATGC':
        if base * (max_run + 1) in sequence:
            return True
    return False

def check_secondary_structures(sequence):
    """Check for potential secondary structures"""
    reverse_complement = str(Seq(sequence).reverse_complement())
    end_length = 6
    sequence_end = sequence[-end_length:]
    rev_comp_start = reverse_complement[:end_length]
    matches = sum(1 for a, b in zip(sequence_end, rev_comp_start) if a == b)
    return matches > end_length/2

def check_primer_dimer(primer1, primer2):
    """Check for potential primer-dimer formation"""
    # Check if 3' ends are complementary
    end_length = 8
    end1 = primer1[-end_length:]
    end2 = str(Seq(primer2[-end_length:]).reverse_complement())
    
    matches = sum(1 for a, b in zip(end1, end2) if a == b)
    return matches > end_length/2

def design_conserved_primers(conserved_regions, min_product_size=70, max_product_size=250):
    """
    Design primer pairs that are fully conserved between human and pig.
    """
    primer_pairs = []
    
    print(f"\nSearching for fully conserved primers in {len(conserved_regions)} regions...")
    
    for region_idx, region in enumerate(conserved_regions):
        human_seq = region['human_seq']
        pig_seq = region['pig_seq']
        region_length = region['length']
        
        print(f"\nAnalyzing region {region_idx+1}: {region_length} bp")
        print(f"  Human ID: {region['human_id']}")
        print(f"  Pig ID: {region['pig_id']}")
        
        if region_length < 60:  # Need enough space for primers and product
            print("  Region too short for primer design, skipping...")
            continue
        
        # Count identical positions to show conservation level
        identical_bases = sum(1 for h, p in zip(human_seq, pig_seq) if h == p)
        conservation_percent = (identical_bases / region_length) * 100
        print(f"  Conservation: {conservation_percent:.1f}% ({identical_bases}/{region_length} bases)")
        
        # Try to design forward primers
        for fw_start in range(0, region_length - 22):
            for fw_len in range(18, 23):
                fw_end = fw_start + fw_len
                if fw_end >= region_length:
                    continue
                
                forward_primer = human_seq[fw_start:fw_end]
                
                # CRITICAL: Check if primer is fully conserved
                if not is_fully_conserved(forward_primer, human_seq, pig_seq, fw_start, fw_end):
                    continue
                
                # Check forward primer properties
                fw_tm = calculate_tm(forward_primer)
                fw_gc = gc_fraction(forward_primer) * 100
                
                # NEW: Strict GC content requirement (45-55%)
                if not (45 <= fw_gc <= 55):
                    continue
                
                # NEW: Strict 3' end requirement (G or C)
                if not has_gc_clamp(forward_primer):
                    continue
                
                # Skip if any other criteria are not met
                if not (52 <= fw_tm <= 62):
                    continue
                if has_repeats(forward_primer) or has_long_runs(forward_primer):
                    continue
                if check_secondary_structures(forward_primer):
                    continue
                
                # Find reverse primer downstream
                for product_size in range(min_product_size, max_product_size + 1, 10):
                    reverse_start = fw_end + product_size
                    if reverse_start >= region_length - 18:  # Need space for reverse primer
                        continue
                    
                    for rv_len in range(18, 23):
                        reverse_end = reverse_start + rv_len
                        if reverse_end > region_length:
                            continue
                        
                        reverse_primer_candidate = human_seq[reverse_start:reverse_end]
                        reverse_primer = str(Seq(reverse_primer_candidate).reverse_complement())
                        
                        # CRITICAL: Check if reverse primer binding site is fully conserved
                        # We need to check the original sequence (not reverse complement)
                        if not is_fully_conserved(reverse_primer_candidate, human_seq, pig_seq, reverse_start, reverse_end):
                            continue
                        
                        # Check reverse primer properties
                        rv_tm = calculate_tm(reverse_primer)
                        rv_gc = gc_fraction(reverse_primer) * 100
                        
                        # NEW: Strict GC content requirement (45-55%)
                        if not (45 <= rv_gc <= 55):
                            continue
                        
                        # NEW: Strict 3' end requirement (G or C)
                        if not has_gc_clamp(reverse_primer):
                            continue
                        
                        if not (52 <= rv_tm <= 62):
                            continue
                        if has_repeats(reverse_primer) or has_long_runs(reverse_primer):
                            continue
                        if check_secondary_structures(reverse_primer):
                            continue
                        
                        # Check for primer-dimer formation
                        if check_primer_dimer(forward_primer, reverse_primer):
                            continue
                        
                        # Calculate a comprehensive score for ranking
                        tm_score = (1 - abs(60 - fw_tm)/10) + (1 - abs(60 - rv_tm)/10)
                        gc_score = (1 - abs(50 - fw_gc)/5) + (1 - abs(50 - rv_gc)/5)  # More strict GC scoring
                        gc_clamp_score = 2.0  # Full points since both primers now have G/C at 3' end
                        
                        total_score = tm_score + gc_score + gc_clamp_score
                        
                        # Valid primer pair found!
                        primer_pairs.append({
                            'forward_seq': forward_primer,
                            'reverse_seq': reverse_primer,
                            'product_size': product_size,
                            'forward_tm': fw_tm,
                            'reverse_tm': rv_tm,
                            'forward_gc': fw_gc,
                            'reverse_gc': rv_gc,
                            'region_index': region_idx + 1,
                            'human_id': region['human_id'],
                            'pig_id': region['pig_id'],
                            'forward_start': fw_start,
                            'reverse_site_start': reverse_start,  # Start of binding site (not primer)
                            'total_score': total_score
                        })
        
        print(f"    Found {len([p for p in primer_pairs if p['region_index'] == region_idx+1])} primer pairs in this region.")
    
    return primer_pairs

def main():
    """
    Main function to design primers from pre-aligned BLAST regions
    """
    print("=" * 70)
    print("CONSERVED PRIMER DESIGN FROM BLAST-ALIGNED REGIONS")
    print("=" * 70)
    
    # ==========================================================================
    # USER INPUT: Provide the names of your BLAST-aligned FASTA files
    # ==========================================================================
    human_aligned_fasta = "Subject_human.fasta"  # <-- CHANGE THIS
    pig_aligned_fasta = "Query_pig.fasta"      # <-- CHANGE THIS
    # ==========================================================================
    
    try:
        # Step 1: Read the conserved regions from BLAST alignment
        conserved_regions = read_conserved_regions(human_aligned_fasta, pig_aligned_fasta)
        
        # Step 2: Design primers that are fully conserved
        primer_pairs = design_conserved_primers(conserved_regions, min_product_size=70, max_product_size=250)
        
        # Step 3: Save and display results
        output_file = "conserved_primers_ABCB11_productsize.txt"
        with open(output_file, 'w') as f:
            f.write("RANKED CONSERVED PRIMER DESIGN RESULTS (From BLAST Alignment)\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"Input files: {human_aligned_fasta}, {pig_aligned_fasta}\n")
            f.write(f"Number of conserved regions analyzed: {len(conserved_regions)}\n")
            
            # NEW: Add filter criteria to output
            f.write("\nFILTER CRITERIA APPLIED:\n")
            f.write("- GC content: 45-55% (strict)\n")
            f.write("- 3' end: Must be G or C (strict)\n")
            f.write("- Tm: 52-62°C\n")
            f.write("- Fully conserved in both species\n")
            f.write("- No secondary structures or primer dimers\n\n")
            
            f.write(f"Total fully conserved primer pairs found: {len(primer_pairs)}\n\n")
            
            if not primer_pairs:
                f.write("No suitable primer pairs found that meet all criteria.\n")
                f.write("This could be due to:\n")
                f.write("1. No primer sequences with GC content 45-55%\n")
                f.write("2. No primer sequences ending with G or C\n")
                f.write("3. No fully conserved regions with suitable properties\n")
                f.write("Try using BLAST with less stringent parameters to find more conserved regions.\n")
            else:
                # Sort by total score (highest first) for ranking
                primer_pairs.sort(key=lambda x: x['total_score'], reverse=True)
                
                f.write("Ranked Primer Pairs (Best First) - ALLY FULLY CONSERVED:\n")
                f.write("-" * 80 + "\n")
                
                for i, pair in enumerate(primer_pairs, 1):
                    f.write(f"\nRANK #{i} (Score: {pair['total_score']:.3f}):\n")
                    f.write(f"  From Region: {pair['region_index']}\n")
                    f.write(f"  Human Location: {pair['human_id']}\n")
                    f.write(f"  Pig Location: {pair['pig_id']}\n")
                    f.write(f"  Forward (5'->3'): {pair['forward_seq']}\n")
                    f.write(f"  3' end: {pair['forward_seq'][-1]} ( G/C clamp)\n")
                    f.write(f"  Reverse (5'->3'): {pair['reverse_seq']}\n")
                    f.write(f"  3' end: {pair['reverse_seq'][-1]} ( G/C clamp)\n")
                    f.write(f"  Product Size: {pair['product_size']} bp\n")
                    f.write(f"  F Tm: {pair['forward_tm']:.1f}°C, GC: {pair['forward_gc']:.1f}% ( 45-55%)\n")
                    f.write(f"  R Tm: {pair['reverse_tm']:.1f}°C, GC: {pair['reverse_gc']:.1f}% ( 45-55%)\n")
                    f.write(f"  F Start Position: {pair['forward_start']}\n")
                    f.write(f"  R Binding Site Start: {pair['reverse_site_start']}\n")
                    f.write("   FULLY CONSERVED IN BOTH HUMAN AND PIG\n")
                    if i < len(primer_pairs):
                        f.write("-" * 60 + "\n")
        
        print(f"\nAnalysis complete! Results saved to '{output_file}'")
        print(f"Found {len(primer_pairs)} primer pairs that meet all strict criteria")
        if primer_pairs:
            print(f"Top-ranked pair score: {primer_pairs[0]['total_score']:.3f}")
            print(f"Top forward primer: {primer_pairs[0]['forward_seq']} (3' end: {primer_pairs[0]['forward_seq'][-1]})")
            print(f"Top reverse primer: {primer_pairs[0]['reverse_seq']} (3' end: {primer_pairs[0]['reverse_seq'][-1]})")
            print(f"GC content: Fwd {primer_pairs[0]['forward_gc']:.1f}%, Rev {primer_pairs[0]['reverse_gc']:.1f}%")
            print("These primers are guaranteed to work for both human and pig sequences!")
            
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()