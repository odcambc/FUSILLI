import pysam

from cigar import Cigar

input_bam = snakemake.input[0]

domain_dict = snakemake.params["domain_dict"]
domain_list = snakemake.params["domain_list"]


def parse_bam_single(input_file):
    # initialize the single counts lists
    single_counts_dict = {}
    for domain in domain_list:
        single_counts_dict[domain] = 0

    # Now make the dict for fusions, indexed by first domain, second domain, then counts of fusion positions

    bamfile = pysam.AlignmentFile(input_file, "rb")

    size = os.path.getsize(bamfile)

    for bam_line in bamfile:
        fusion = 0
        line = bam_line.tostring().split("\t")
        for field in line:
            if field.startswith("SA:"):
                fusion = 1

        if not fusion:
            if line[2] not in single_counts_dict.keys():
                single_counts_dict[line[2]] = 0

            single_counts_dict[line[2]] += 1

    return


def parse_bam_fusion(input_file):
    """Parse the fusion reads from the filtered, sorted bam file.
    All of the reads in the input file have an SA tag, and so should have
    supplementary alignments. We don't need to check for its presence.

    The function will parse the bam file and identify, for each read,
    the domain of the first alignment and the domain of the supplementary alignment.
    It will use the starting and end positions of the alignments using the cigar strings.
    It will then check if the read covers the breakpoint by checking if the end position of
    the alignment, for the full-length domain, extends to its end.

    The output is a csv file with the following columns:
    domain1: The domain of the first alignment
    domain2: The domain of the supplementary alignment
    pos1: The starting position of the first alignment
    end_pos1: The end position of the first alignment
    pos2: The starting position of the supplementary alignment
    end_pos2: The end position of the supplementary alignment
    full1: A boolean indicating if the read covers the full extent of the first domain
    full2: A boolean indicating if the read covers the full extent of the second domain
    count: The number of reads with the given alignment positions and domains.
    """

    fusion_counts_dict = {}

    match = 0

    bamfile = pysam.AlignmentFile(input_file, "rb")

    for bam_line in bamfile:
        line = bam_line.tostring().split("\t")
        for field in line:
            if field.startswith("SA:"):
                supplemental_field = field
                match = 1
                break

        if match:
            # sequence = line[9]
            # The primary alignment is in the main fields, as standard SAM format.
            chrom = line[2]
            pos = line[3]
            cigar = line[5]

            # The supplementary alignment is in the SA field in the tags.
            # Example:
            # SA:Z:YWHAE,680,+,35M116S,1,0;
            # Split by colon, then get the 3rd field, which is the supplementary alignment.
            # Split by comma, then get the 1st, 2nd, 4th fields, which are chrom, pos, cigar of the supplementary alignment.
            supplemental_fields = supplemental_field.split(":")[2].split(",")

            chrom_supplemental = supplemental_fields[0]
            pos_supplemental = supplemental_fields[1]
            cigar_supplemental = supplemental_fields[3]

            cigar_len = get_cigar_length(cigar)

            end_pos = int(pos) + cigar_len

            # For the supplementary alignment, parse the cigar string.

            cigar_len = get_cigar_length(cigar_supplemental)

            end_pos_supplemental = int(pos_supplemental) + cigar_len

            # Add the count to the dict

            # Check whether the domains are full length.

            full = 0
            full2 = 0

            if end_pos == domain_dict[chrom]:
                full = 1

            if end_pos_supplemental == domain_dict[chrom_supplemental]:
                full2 = 1

            if (
                chrom,
                chrom_supplemental,
                pos,
                end_pos,
                pos_supplemental,
                end_pos_supplemental,
                full,
                full2,
            ) not in fusion_counts_dict:
                fusion_counts_dict[
                    (
                        chrom,
                        chrom_supplemental,
                        pos,
                        end_pos,
                        pos_supplemental,
                        end_pos_supplemental,
                        full,
                        full2,
                    )
                ] = 1
            else:
                fusion_counts_dict[
                    (
                        chrom,
                        chrom_supplemental,
                        pos,
                        end_pos,
                        pos_supplemental,
                        end_pos_supplemental,
                        full,
                        full2,
                    )
                ] += 1

            # Reset the match flag
            match = 0

    # Write the dict.

    write_bam_fusion(fusion_counts_dict, snakemake.output[0])

    return


def get_cigar_length(cigar):
    cigar_list = list(Cigar(cigar).items())

    cigar_len = 0

    for cigar_tuple in cigar_list:
        match cigar_tuple[1]:
            case "M":
                cigar_len += cigar_tuple[0]
            case "X":
                cigar_len += cigar_tuple[0]
            case "D":
                cigar_len += cigar_tuple[0]
            case "N":
                cigar_len += cigar_tuple[0]
            case "I":
                cigar_len -= cigar_tuple[0]
            case "S":
                pass
            case "H":
                pass
            case "P":
                pass
            case _:
                raise ValueError(f"Unexpected cigar operation: {cigar_tuple[1]}")

    cigar_len = cigar_len - 1

    return cigar_len


def write_bam_fusion(fusion_counts_dict, output_csv):
    # Keys in the dict are tuples of the form
    # (chrom, chrom_supplemental, pos, end_pos, pos_supplemental, end_pos_supplemental, cigar, cigar_supplemental, full, full2)
    # The values are the counts of the fusions.

    with open(output_csv, "w", encoding="utf-8") as f:
        f.write("domain1,domain2,pos1,end_pos1,pos2,end_pos2,full,full2,count\n")
        for key, value in fusion_counts_dict.items():
            f.write(
                f"{key[0]},{key[1]},{key[2]},{key[3]},{key[4]},{key[5]},{key[6]},{key[7]},{value}\n"
            )
    return


print("Parsing BAM file: ", input_bam)
parse_bam_fusion(input_bam)
# parse_bam_single(input_bam)
print("Done parsing BAM file: ", input_bam)
