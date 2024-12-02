from Bio import SeqIO
from os import listdir
import argparse

def genome_statistics(fasta_dir, fasta_suffix):
    '''
    INPUT: A directory containing FASTA files and the suffix of the FASTA files (e.g., .fas)
    OUTPUT: A CSV file with the genome size, number of contigs, and N50 for each FASTA file
    '''
    fasta_stats = {}

    for fasta_file in listdir(fasta_dir):
        if fasta_file.endswith(fasta_suffix):
            contig_lengths = []
            for seq_record in SeqIO.parse(f"{fasta_dir}/{fasta_file}", 'fasta'):
                contig_lengths.append(len(seq_record.seq))

            if contig_lengths:
                contig_lengths.sort(reverse=True)
                genome_size = sum(contig_lengths)
                num_contigs = len(contig_lengths)

                # Calculate N50
                cumulative_length = 0
                n50 = 0
                for length in contig_lengths:
                    cumulative_length += length
                    if cumulative_length >= genome_size / 2:
                        n50 = length
                        break

                fasta_stats[fasta_file] = {
                    'genome_size': genome_size,
                    'num_contigs': num_contigs,
                    'n50': n50
                }

    with open('fasta_stats.csv', 'w') as f:
        f.write("File,Genome_Size,Num_Contigs,N50\n")
        for fasta_file, stats in fasta_stats.items():
            f.write(f"{fasta_file},{stats['genome_size']},{stats['num_contigs']},{stats['n50']}\n")

def main(args=None):
    parser = argparse.ArgumentParser(description='Calculate genome statistics for FASTA files.')
    parser.add_argument('--input', '-i', help='Input directory containing FASTA files.')
    parser.add_argument('--suffix', '-s', help='Suffix of the FASTA files (e.g., .fasta, .fas).')
    args = parser.parse_args(args)

    genome_statistics(args.input, args.suffix)

if __name__ == "__main__":
    main()
