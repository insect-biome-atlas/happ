#!/usr/bin/env python

from argparse import ArgumentParser


def main(args):
    header = "\t".join(["ASV"] + args.ranks)
    with open(args.infile, "r") as fhin, open(args.outfile, "w") as fhout:
        for i, line in enumerate(fhin):
            if i == 0:
                fhout.write(f"{header}\n")
                continue
            name, lwr, fract, alwr, afract, taxopath = line.rstrip().split("\t")
            lwr = float(lwr)
            fract = float(fract)
            alwr = float(alwr)
            afract = float(afract)
            taxnames = taxopath.split(";")
            taxnames += [f"unclassified.{taxnames[-1]}"] * (
                    len(args.ranks) - len(taxnames)
            )
            line_out = "\t".join([name] + taxnames)
            fhout.write(f"{line_out}\n")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("infile", type=str,
                        help="Gappa best-hit results")
    parser.add_argument("outfile", type=str,
                        help="Output file name")
    parser.add_argument("--ranks", nargs="+", default=["order","family","genus","species"],
                        help="Ranks to write to output")
    args = parser.parse_args()
    main(args)