import argparse
from Bio import SeqIO

def fasta_parsing(fasta, fasta_dict, tag):
    fasta_file = SeqIO.parse(fasta, "fasta")
    for el in fasta_file:
        name, sequence = el.id, el.seq
        if tag == "nucl":
            fasta_dict[name] = sequence
        elif tag == "prot":
            fasta_dict[name.split(".p")[0]] = {"protein_ID": name, "sequence": sequence,
                                               "length": len(sequence)}


def expression_table_parsing(exp_table, header, exp_dict):
    header.extend(exp_table.readline().strip().split("\t"))
    print("Header of the table with exp.values: {header}".format(header='\t'.join(header)))
    for line in exp_table:
        description = line.strip().split("\t")
        gene, values = description[0], description[1:]
        exp_dict[gene] = {sample: values[header.index(sample) - 1] for sample in header[1:]}


def map_parsing(map_dict, gene_trans_map):
    for line in gene_trans_map:
        description = line.strip().split('\t')
        gene, trans = description[0], description[1]
        if gene not in map_dict.keys():
            map_dict[gene] = []
            map_dict[gene].append(trans)
        else:
            map_dict[gene].append(trans)


def filters(protein_dict, averaged_exp_dict, map_dict, min_len, wanted_dict, out):
    #   FIRST FILTER: noticeable (AVERAGED) expression level
    noticeable_expression_list = []
    noticeable_expression_list.extend(averaged_exp_dict.keys())
    #   SECOND FILTER: protein-coding and protein`s length >= min_len
    with open("{out}.filters_log".format(out=out), 'a') as output:
        output.write("Gene_ID\tLengths\tChose\tTranscript_ID\n")
        for gene in noticeable_expression_list:
            protein_sizes_list, trans_ids_list = [], []
            for trans in map_dict[gene]:
                if trans in protein_dict.keys():
                    protein_sizes_list.append(protein_dict[trans]["length"])
                    trans_ids_list.append(trans)

            if len(protein_sizes_list) != 0 and max(protein_sizes_list) >= min_len:
                for trans in trans_ids_list:
                    if protein_dict[trans]["length"] == max(protein_sizes_list):
                        output.write("{gene}\t{len}\t{max}\t{id}\n".format(
                            gene=gene, len=";".join([str(el) for el in protein_sizes_list]),
                            max=protein_dict[trans]["length"], id=trans))
                        wanted_dict[gene] = {"transcript_ID": trans, "protein_ID": protein_dict[trans]["protein_ID"]}
                        break

    print("After first filter (noticeable averaged expression level) survived: {first}\n"
          "After second filter (protein-coding sequence and protein`s length >= {min} aa) survived: {second}".format(
            first=len(noticeable_expression_list), min=min_len, second=len(wanted_dict.keys())))


def output_writing(out, nucl_dict, prot_dict, unaveraged_header, unaveraged_exp_dict,
                   averaged_header, averaged_exp_dict, wanted_dict):
    wanted_genes, wanted_transcripts, wanted_proteins = [], [], []
    with open("{out}.gene_map.tsv".format(out=out), 'a') as all_wanted_output:
        all_wanted_output.write("Genes\tTranscript_IDs\tProtein_IDs\n")
        for gene, values in wanted_dict.items():
            wanted_genes.append(gene)
            wanted_transcripts.append(values["transcript_ID"])
            wanted_proteins.append(values["protein_ID"])
            all_wanted_output.write("{gene}\t{transcript}\t{protein}\n".format(gene=gene,
                                                                               transcript=values["transcript_ID"],
                                                                               protein=values["protein_ID"]))

    with open("{out}.genes_level.after_filters.nucl.fasta".format(out=out), 'a') as output_nucl:
        for contig, sequence in nucl_dict.items():
            if contig in wanted_transcripts:
                output_nucl.write(">{name}\n{seq}\n".format(name=contig, seq=sequence))

    with open("{out}.genes_level.after_filters.prot.fasta".format(out=out), 'a') as output_prot:
        for protein, values in prot_dict.items():
            if values["protein_ID"] in wanted_proteins:
                output_prot.write(">{name}\n{seq}\n".format(name=values["protein_ID"], seq=values["sequence"]))

    with open("{out}.genes_level.after_filters.unaveraged_TPMs.tsv".format(out=out), 'a') as output_unaveraged:
        output_unaveraged.write("{header}\n".format(header="\t".join(unaveraged_header)))
        for gene, values in unaveraged_exp_dict.items():
            if gene in wanted_genes:
                ordered_values = [values[sample] for sample in unaveraged_header[1:]]
                output_unaveraged.write("{gene}\t{ordered_values}\n".format(gene=gene, ordered_values="\t".join(
                    ordered_values
                )))

    with open("{out}.genes_level.after_filters.averaged_TPMs.tsv".format(out=out), 'a') as output_averaged:
        output_averaged.write("{header}\n".format(header="\t".join(averaged_header)))
        for gene, values in averaged_exp_dict.items():
            if gene in wanted_genes:
                ordered_values = [values[sample] for sample in averaged_header[1:]]
                output_averaged.write("{gene}\t{ordered_values}\n".format(gene=gene, ordered_values="\t".join(
                    ordered_values)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--nucl_fasta', type=argparse.FileType('r'), required=True,
                        help="Fasta file with results of clusterization performed with CDHIT-est software")
    parser.add_argument('--prot_fasta', type=argparse.FileType('r'), required=True,
                        help="Fasta file with aminoacid sequences predicted by TransDecoder software")
    parser.add_argument('--gene_trans_map', type=argparse.FileType('r'), required=True,
                        help="File with 'gene_trans_map' description created by Trinity")
    parser.add_argument('--unaveraged_exp', type=argparse.FileType('r'), required=True,
                        help="Table with unaveraged expression levels of selected GENES")
    parser.add_argument('--averaged_exp', type=argparse.FileType('r'), required=True,
                        help="Table with averaged expression levels of selected GENES "
                             "(for instance, genes with averaged between replicates TPM > 1)")
    parser.add_argument('--min_len', type=int, required=True,
                        help="Minimal protein length in amino acids (for instance, 75 or 100)")
    parser.add_argument('--out', type=str, required=True)
    args = parser.parse_args()
    
    nucl_dict, prot_dict = {}, {}
    unaveraged_header, unaveraged_exp_dict = [], {}
    averaged_header, averaged_exp_dict = [], {}
    map_dict, wanted_dict = {}, {}

    print("##### Input files parsing #####")
    fasta_parsing(args.nucl_fasta, nucl_dict, "nucl")
    fasta_parsing(args.prot_fasta, prot_dict, "prot")
    expression_table_parsing(args.unaveraged_exp, unaveraged_header, unaveraged_exp_dict)
    expression_table_parsing(args.averaged_exp, averaged_header, averaged_exp_dict)
    map_parsing(map_dict, args.gene_trans_map)
    print("##### Applying filters #####")
    filters(prot_dict, averaged_exp_dict, map_dict, args.min_len, wanted_dict, args.out)
    print("##### Output files creation #####")
    output_writing(args.out, nucl_dict, prot_dict, unaveraged_header, unaveraged_exp_dict,
                   averaged_header, averaged_exp_dict, wanted_dict)

