#!/usr/bin/env python

from argparse import ArgumentParser
from subprocess import Popen, check_output, PIPE
from pathlib import Path


def check_lines(f):
    command = ["wc", "-l", f]
    l = check_output(command)
    lines = int(l.decode().lstrip(" ").split(" ")[0])
    return lines


def handle_empty(outdir, counts):
    with open(f"{outdir}/asv_seqs.opti_mcc.list", "w") as fhout, open(
        counts, "r"
    ) as fhin:
        fhout.write("No results\n")
        for i, line in enumerate(fhin):
            if i == 0:
                continue
            asv = line.rstrip().split("\t")[0]
            fhout.write(f"{asv}\n")
    Path(f"{outdir}/asv_seqs.opti_mcc.sensspec").touch()
    Path(f"{outdir}/asv_seqs.opti_mcc.steps").touch()
    return


def main(args):
    lines = check_lines(args.dist)
    if lines == 0:
        handle_empty(args.outdir, args.counts)
        return
    command = f"""
        mothur "#set.dir(output={args.outdir});cluster(column={args.dist}, count={args.counts},method=opti, delta={args.delta}, cutoff={args.cutoff}, initialize={args.initialize},precision={args.precision})"
        """
    print(f"Running {command}")
    p = Popen(command, stderr=PIPE, stdout=PIPE, shell=True)
    p.wait()


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("dist", type=str, help="Input distance file")
    parser.add_argument("counts", type=str, help="Counts file")
    parser.add_argument("--delta", type=float, help="Delta param")
    parser.add_argument("--cutoff", type=float, help="Cutoff param")
    parser.add_argument("--initialize", type=str, help="Initialize param")
    parser.add_argument("--precision", type=int, help="Precision param")
    parser.add_argument("--outdir", type=str, help="Output directory")
    parser.add_argument("--log", type=str, help="Mothur logfile")
    args = parser.parse_args()
    main(args)
