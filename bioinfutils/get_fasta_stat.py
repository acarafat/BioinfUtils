from Bio import SeqIO
import argparse

def fasta_statistics(fasta_file):
    '''
    INPUT: Path to a FASTA file
    OUTPUT: Print output (contig, contig size) on the screen 
    '''

    for seq_record in SeqIO.parse(fasta_file, 'fasta'):
            print( f"{seq_record.id} has a length of: {len(seq_record.seq)} bp")

    pass

def main(args=None):
    parser = argparse.ArgumentParser(description='Print statistic of contigs on a fasta file')
    parser.add_argument('--input', '-i', help='Path to FASTA files')
    args = parser.parse_args(args)

    fasta_statistics(args.input)

if __name__ == "__main__":
    main()
