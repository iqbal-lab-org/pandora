import argparse

def sam_to_fasta(sam_fh, fasta_fh):
    for line in sam_fh:
        is_header = line.startswith("@")
        if is_header:
            continue
        line_split = line.split("\t")
        is_mapped = line_split[1] == "0"
        if is_mapped:
            print(f">{line_split[0]}~~~{line_split[2]}~~~{line_split[5]}", file=fasta_fh)
            print(line_split[9], file=fasta_fh)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sam", action="store", type=str, required=True,
                        help='Pandora sam file')
    parser.add_argument("--fasta", action="store", type=str, required=True,
                        help='Ouput fasta file')
    args = parser.parse_args()

    with open(args.sam) as sam_fh, open(args.fasta, "w") as fasta_fh:
        sam_to_fasta(sam_fh, fasta_fh)

main()