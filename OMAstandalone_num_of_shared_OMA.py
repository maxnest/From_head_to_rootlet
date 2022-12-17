import argparse

def matrix_parsing(matrix, matrix_dict, species_list):
    species_list.extend(matrix.readline().strip().split("\t")[1:])
    print("Number of species: {num}\nList of species: {species}\n".format(
        num=len(species_list), species="\t".join(species_list)))
    for species in species_list:
        matrix_dict[species] = []

    for line in matrix:
        description = line.strip().split("\t")
        oma, values = description[0], description[1:]
        for i, el in enumerate(values):
            if el == "1":
                matrix_dict[species_list[i]].append(oma)


def species_comparison(species_list, results_dict):
    for species in species_list:
        for other_species in species_list:
            pair_tag = "{sp}|{other}".format(sp=species, other=other_species)
            results_dict[pair_tag] = set.intersection(*[set(matrix_dict[species]), set(matrix_dict[other_species])])


def output_writing(out, species_list, result_dict):
    with open("{out}.num_shared_OMAs.tsv".format(out=out), 'a') as output_file:
        output_file.write("Species\Species\t{species}\n".format(species="\t".join(species_list)))
        for species in species_list:
            shared_oma = []
            for other_species in species_list:
                shared_oma.append(str(len(result_dict["{sp}|{other}".format(sp=species, other=other_species)])))
            output_file.write("{species}\t{shared_oma}\n".format(species=species, shared_oma="\t".join(shared_oma)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--matrix', type=argparse.FileType('r'), required=True,
                        help="Presence/absence matrix created by OMAstandalone (PhyleticProfileOMAGroups.txt) "
                             "without first rows with description")
    parser.add_argument('--out', type=str, required=True)
    args = parser.parse_args()
    
    matrix_dict, species_list, result_dict = {}, [], {}
    matrix_parsing(args.matrix, matrix_dict, species_list)
    species_comparison(species_list, result_dict)
    output_writing(args.out, species_list, result_dict)