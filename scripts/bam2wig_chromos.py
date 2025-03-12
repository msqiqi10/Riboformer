import pysam
import argparse
from collections import defaultdict

def extract_chromosome_name(reference_name):
    """
    Extract chromosome identifier from the reference name.

    Args:
        reference_name (str): Reference name from BAM file (e.g., 'AT1G01010.1').

    Returns:
        str: Chromosome identifier matching those in seq_dict (e.g., '1', '2', '3', 'mt', 'pt').
    """
    # Adjust parsing logic based on your naming conventions
    if reference_name.startswith('Chr'):
        return reference_name.replace('Chr', '')
    elif reference_name.startswith('AT'):
        # Assuming 'AT1G' corresponds to chromosome '1'
        return reference_name[2]
    else:
        return reference_name  # Return as-is if it already matches

def convert_bam_to_wig(bam_file, output_prefix, p_site_offset=14):
    """
    Convert BAM file to forward and reverse WIG files with P-site offset.

    Args:
        bam_file (str): Input BAM file path.
        output_prefix (str): Prefix for output WIG files.
        p_site_offset (int): P-site offset from 5' end of read.
    """
    # Initialize counters for forward and reverse strands per chromosome
    forward_counts = defaultdict(lambda: defaultdict(int))  # {chrom: {pos: count}}
    reverse_counts = defaultdict(lambda: defaultdict(int))  # {chrom: {pos: count}}

    # Open BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")

    # Store all chromosomes
    chromosomes = bam.references
    chr_lengths = dict(zip(bam.references, bam.lengths))

    # Process each read
    for read in bam.fetch():
        # Skip unmapped reads
        if read.is_unmapped:
            continue

        # Get chromosome name and adjust it to match seq_dict keys
        chrom_full = bam.get_reference_name(read.reference_id)
        chrom = extract_chromosome_name(chrom_full)

        # Calculate P-site position
        if read.is_reverse:
            # Adjust for 0-based indexing and exclusive end coordinate
            p_site_pos = read.reference_end - 1 - p_site_offset
            reverse_counts[chrom][p_site_pos] += 1
        else:
            p_site_pos = read.reference_start + p_site_offset
            forward_counts[chrom][p_site_pos] += 1

    bam.close()

    # Write forward strand WIG file
    with open(f"{output_prefix}_f.wig", "w") as f_out:
        f_out.write("track type=wiggle_0 name=forward\n")
        for chrom in chromosomes:
            chrom_adjusted = extract_chromosome_name(chrom)
            chrom_length = chr_lengths[chrom]
            counts = forward_counts.get(chrom_adjusted, {})
            if counts:
                f_out.write(f"fixedStep chrom={chrom_adjusted} start=1 step=1\n")
                for pos in range(chrom_length):
                    count = counts.get(pos, 0)
                    f_out.write(f"{count}\n")
            else:
                # Write zeros if no counts for this chromosome
                f_out.write(f"fixedStep chrom={chrom_adjusted} start=1 step=1\n")
                for _ in range(chrom_length):
                    f_out.write("0\n")

    # Write reverse strand WIG file
    with open(f"{output_prefix}_r.wig", "w") as r_out:
        r_out.write("track type=wiggle_0 name=reverse\n")
        for chrom in chromosomes:
            chrom_adjusted = extract_chromosome_name(chrom)
            chrom_length = chr_lengths[chrom]
            counts = reverse_counts.get(chrom_adjusted, {})
            if counts:
                r_out.write(f"fixedStep chrom={chrom_adjusted} start=1 step=1\n")
                for pos in range(chrom_length):
                    count = counts.get(pos, 0)
                    r_out.write(f"{count}\n")
            else:
                # Write zeros if no counts for this chromosome
                r_out.write(f"fixedStep chrom={chrom_adjusted} start=1 step=1\n")
                for _ in range(chrom_length):
                    r_out.write("0\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert BAM to strand-specific WIG files')
    parser.add_argument('-i', '--input', required=True, help='Input BAM file')
    parser.add_argument('-o', '--output_prefix', required=True, help='Output prefix for WIG files')
    parser.add_argument('-p', '--p_site_offset', type=int, default=14,
                        help='P-site offset from 5\' end of read (default: 14)')

    args = parser.parse_args()

    convert_bam_to_wig(args.input, args.output_prefix, args.p_site_offset)