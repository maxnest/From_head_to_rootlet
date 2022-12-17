import argparse
from Bio import SeqIO

def fasta_parsing(gene_map, fasta_seqs):
    for fasta in fasta_seqs:
        gene_map["g{num}".format(num=len(gene_map.keys()) + 1)] = {"name": fasta.id, "seq": fasta.seq}


def output_files_creating(gene_map, output):
    with open("{output}.renamed_for_secretomeP.fasta".format(output=output), 'a') as fasta_output:
        for gene, values in gene_map.items():
            fasta_output.write(">{gene}\n{seq}\n".format(gene=gene, seq=values["seq"]))

    with open("{output}.old_2_new_seq_name_map.tsv".format(output=output), 'a') as map_output:
        map_output.write("Old_name\tNew_name\n")
        for gene, values in gene_map.items():
            map_output.write("{old}\t{new}\n".format(old=values["name"], new=gene))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=argparse.FileType('r'), required=True)
    parser.add_argument('--output', type=str)
    args = parser.parse_args()
    
    gene_map = {}
    fasta_seqs = SeqIO.parse(args.fasta, "fasta")
    fasta_parsing(gene_map, fasta_seqs)
    output_files_creating(gene_map, args.output)