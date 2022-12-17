import argparse
from Bio import SeqIO


def tmhmm_parsing(tmhmm, with_tmhelix, without_tmhelix):
    tmhmm_dict = {}
    for line in tmhmm:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            ID, prediction = description[0], description[2]
            if ID not in tmhmm_dict.keys():
                tmhmm_dict[ID] = []
                tmhmm_dict[ID].append(prediction)
            elif ID in tmhmm_dict.keys():
                tmhmm_dict[ID].append(prediction)

    for ID, predictions in tmhmm_dict.items():
        if "TMhelix" in predictions:
            with_tmhelix.append(ID)
        else:
            without_tmhelix.append(ID)

    print("TMhelix was discovered in {with_tmhelix} sequences, "
          "while there are no such structures in {without_tmhelix}".format(
           with_tmhelix=len(with_tmhelix), without_tmhelix=len(without_tmhelix)))


def fasta_parsing(fasta, output, without_tmhelix):
    fasta_seqs = SeqIO.parse(fasta, "fasta")
    with open("{output}.fasta".format(output=output), 'a') as output_file:
        for fasta in fasta_seqs:
            name, sequence = fasta.id, fasta.seq
            if name in without_tmhelix:
                output_file.write(">{name}\n{seq}\n".format(name=name, seq=sequence))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=argparse.FileType('r'), required=True)
    parser.add_argument('--tmhmm', type=argparse.FileType('r'), required=True)
    parser.add_argument('--output', type=str)
    args = parser.parse_args()
    
    with_tmhelix, without_tmhelix = [], []
    print("***** TMHMM-output parsing *****")
    tmhmm_parsing(args.tmhmm, with_tmhelix, without_tmhelix)
    print("***** Fasta-file parsing *****")
    fasta_parsing(args.fasta, args.output, without_tmhelix)
