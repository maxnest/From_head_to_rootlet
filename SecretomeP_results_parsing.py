import argparse
from Bio import SeqIO

def fasta_parsing(fasta_seqs, fasta_dict):
    for fasta in fasta_seqs:
        fasta_dict[fasta.id] = fasta.seq


def old_2_new_name_parsing(old_2_new, names_dict):
    header = old_2_new.readline()
    for line in old_2_new:
        description = line.strip().split("\t")
        old, new = description[0], description[1]
        names_dict[old] = new


def secretomep_parsing(table, secretomep_dict):
    for line in table:
        if not line.startswith("#"):
            description = [el for el in line.strip().split(" ") if len(el) != 0]
            name, nn, warning = description[0], float(description[1].split("\t")[1]), description[-1]
            secretomep_dict[name] = {"nn": nn, "warning": warning}


def classification(secretomep_dict, nonclassic):
    for name, values in secretomep_dict.items():
        if values["nn"] >= 0.9 and values["warning"] == "-":
            nonclassic.append(name)


def output_writing(output, fasta_dict, names_dict, nonclassic):
    with open("{output}.fasta".format(output=output), 'a') as output_fasta:
        for id, seq in fasta_dict.items():
            if names_dict[id] in nonclassic:
                output_fasta.write(">{id}\n{seq}\n".format(id=id, seq=seq))

    with open("{output}.new_IDs_of_all_nonclassic_exsec.txt".format(output=output), 'a') as ids:
        for id in nonclassic:
            ids.write("{id}\n".format(id=id))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--table', type=argparse.FileType('r'), required=True,
                        help="Summary table with all results of classification")
    parser.add_argument('--old_2_new', type=argparse.FileType('r'), required=True,
                        help="Table with old and new names of contigs")
    parser.add_argument('--fasta', type=argparse.FileType('r'), required=True,
                        help="Initial fasta file")
    parser.add_argument('--output', type=str)
    args = parser.parse_args()
    
    fasta_dict, names_dict, secretomep_dict, nonclassic = {}, {}, {}, []
    fasta_seqs = SeqIO.parse(args.fasta, "fasta")
    fasta_parsing(fasta_seqs, fasta_dict)
    old_2_new_name_parsing(args.old_2_new, names_dict)
    secretomep_parsing(args.table, secretomep_dict)
    classification(secretomep_dict, nonclassic)
    output_writing(args.output, fasta_dict, names_dict, nonclassic)
