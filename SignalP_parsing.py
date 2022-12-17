import argparse

def signalp_parsing(signalp, threshold, list_more, list_less):
    for line in signalp:
        if not line.startswith("#"):
            description = line.strip().split("\t")
            contig_ID, other_prob = description[0], float(description[3])
            if other_prob >= threshold:
                list_more.append(contig_ID)
            elif other_prob < threshold:
                list_less.append(contig_ID)


def write_output(list_more, list_less, output, threshold):
    with open("{out}.more_than_{threshold}.txt".format(out=output, threshold=threshold), 'a') as more:
        for contig in list_more:
            more.write("{contig}\n".format(contig=contig))

    with open("{out}.less_than_{threshold}.txt".format(out=output, threshold=threshold), 'a') as less:
        for contig in list_less:
            less.write("{contig}\n".format(contig=contig))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--summary', type=argparse.FileType('r'), required=True,
                        help="SignalP (v5.0b) output file")
    parser.add_argument('--threshold', type=float, required=True.
                        help = "Threshold value for the OTHER-prob metric")
    parser.add_argument('--output', type=str, required=True)
    args = parser.parse_args()
    
    list_with_more, list_with_less = [], []
    print("*** SignalP summary parsing ***")
    signalp_parsing(args.summary, args.threshold, list_with_more, list_with_less)
    print("*** Output files creating ***")
    write_output(list_with_more, list_with_less, args.output, args.threshold)
