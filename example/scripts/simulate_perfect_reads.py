# very simple script to simulate perfect reads
import argparse
from Bio import SeqIO
import random

def get_args():
    parser = argparse.ArgumentParser(description='Simulate perfect reads from a fasta file.')
    parser.add_argument('--fasta_file', type=str, required=True)
    parser.add_argument('--reads_length', type=int, required=True)
    parser.add_argument('--number_of_reads', type=int, required=True)
    args = parser.parse_args()
    return args


def read_fasta_and_get_sequences(fasta_filepath):
    sequences = []

    for record in SeqIO.parse(fasta_filepath, "fasta"):
        sequences.append(str(record.seq))
        sequences.append(str(record.reverse_complement().seq))

    return sequences


def simulate_reads(sequences_to_get_reads_from, reads_length, number_of_reads):
    for read_id in range(number_of_reads):
        sequence = random.choice(sequences_to_get_reads_from)
        start = random.choice(range(0, len(sequence)-reads_length))
        end = start+reads_length
        simulated_read = sequence[start:end]

        assert len(simulated_read) == reads_length
        print(f">simulated_read_{read_id}")
        print(simulated_read)
        print("+")
        print("+"*reads_length)


def main():
    args = get_args()
    sequences_to_get_reads_from = read_fasta_and_get_sequences(args.fasta_file)
    simulate_reads(sequences_to_get_reads_from, args.reads_length, args.number_of_reads)


if __name__ == "__main__":
    main()