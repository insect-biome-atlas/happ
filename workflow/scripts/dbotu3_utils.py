#!/usr/bin/env python

import os
import shutil


def dbotu32tab(sm):
    """
    Creates a mapping table of asv -> cluster
    """
    os.makedirs(sm.params.tmpdir, exist_ok=True)
    with open(sm.input[0], "r") as fhin, open(sm.params.out, "w") as fhout:
        fhout.write("asv\tcluster\n")
        for i, line in enumerate(fhin):
            line = line.rstrip()
            asvs = line.rsplit()[1:]
            for asv in asvs:
                fhout.write(f"{asv}\tcluster{i}\n")
    shutil.move(sm.params.out, sm.output[0])
    os.rmdir(sm.params.tmpdir)


def main(sm):
    toolbox = {"dbotu32tab": dbotu32tab}
    toolbox[sm.rule](sm)


if __name__ == "__main__":
    main(snakemake)
