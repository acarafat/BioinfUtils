from Bio import SeqIO
import argparse


def byList(seqID_list, fasta_target, output):
    '''
    INPUT:
      seqID: list of sequence ID or fasta identifier
      fasta_target: Fasta file containing sequences which need to be filtered
    OUTPUT: Fasta file containing sequences of those ID
    '''
    filtered_entry = []
    for seq_record in SeqIO.parse(fasta_target, 'fasta'):
        if seq_record.id in seqID_list:
            filtered_entry.append(seq_record)
    SeqIO.write(filtered_entry, output, 'fasta')
    pass


def bySize(fasta_target, new_fasta):
    '''
    INPUT: Fasta file containing duplicate sequences of multiple lenght
    ALOGRITHM: Only keep the longest sequences
    Output: Fasta file contiaing unique sequences
    '''
    unique_seqs = {}
    for seq_record in SeqIO.parse(fasta_target, 'fasta'):
        if seq_record.id not in unique_seqs.keys():
            unique_seqs[seq_record.id] = seq_record
        else:
            if len(str(seq_record.seq)) > len(str(unique_seqs[seq_record.id])):
                unique_seqs[seq_record.id] = seq_record
    SeqIO.write(unique_seqs.values(), new_fasta, 'fasta')
    pass

def filter_fasta_by_id(fasta_file, id_file, output_file):
    """Filters a FASTA file by sequence IDs in another file using Biopython.
          Args:
           fasta_file: Path to the FASTA file.
           id_file: Path to the text file containing IDs (one per line).
           output_file: Path to the output FASTA file.
    """

    with open(id_file, 'r') as ids:
        desired_ids = set(line.strip() for line in ids)

    # Use Biopython's SeqIO to iterate through FASTA records
    with open(fasta_file, 'r') as fasta, open(output_file, 'w') as out:
        for record in SeqIO.parse(fasta, 'fasta'):
            if record.id in desired_ids:
            SeqIO.write(record, out, 'fasta')


def main(args=None):
    parser = argparse.ArgumentParser(description='Filter FASTA file by sequence IDs, size, or lis
t')
    parser.add_argument('--fasta','-f',  required = True, help='Path to the FASTA file')
    parser.add_argument('--output', '-o',  required = True, help='Path to the output FASTA file')
    #parser.add_argument('--option', '-n',  required = True, help='Type of filter: 1 for by id in list, 2 for by size, 3 $')
    parser.add_argument('--id', '-i',  required = False, help='Path to the file containing IDs')
    parser.add_argument('--list', '-l', nargs='+',  required =  False, help='List of contig id\'s to be filtered')

    args = parser.parse_args(args)

    # Call the filtering function
    if args.id != None:
        filter_fasta_by_id(args.fasta, args.ids, args.output)
    elif args.list != None:
        byList(args.list, args.fasta, args.output)
    else:
        bySize(args.fasta, args.output)

    print(f"Filtered FASTA written to: {args.output}")

if __name__ == '__main__':
    main()
