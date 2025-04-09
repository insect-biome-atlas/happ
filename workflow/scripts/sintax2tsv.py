#!/usr/bin/env python
import pandas as pd
import re
from argparse import ArgumentParser


def parse_taxonomy(df, r, ranks=["kingdom", "phylum", "class", "order", "family", "genus", "species"], cutoff=0.8):
    """
    Parse the taxonomy from the SINTAX output

    This function should be used with the pandas groupby method, grouping by the
    ASV id.

    Taking a list of ranks and a cutoff value, it will parse the taxonomy string
    and return a dataframe with the parsed taxonomy.
    """
    asv_id = df.index[0]
    tax = df.iloc[0, 0]
    conf_values = r.findall(tax)
    tax_values = [x[2:] for x in r.sub("", tax).split(",")]
    d = {asv_id: {}}
    last_known = ""
    for i, rank in enumerate(ranks):
        try:
            d[asv_id][f"{rank}_conf"] = conf_values[i]
        except IndexError:
            d[asv_id][f"{rank}_conf"] = None
            d[asv_id][rank] = None
            continue
        if float(conf_values[i]) >= cutoff:
            d[asv_id][rank] = tax_values[i]
            d[asv_id][f"{rank}_conf"] = conf_values[i]
            last_known = tax_values[i]
        else:
            d[asv_id][rank] = f"unclassified{'.'*min(len(last_known), 1)}{last_known}"
    return pd.DataFrame(d,).T

def main(infile, outfile, conf_outfile, ranks, cutoff):
    r = re.compile(r"\(([0-9]\.[0-9]+)\)")
    sintax = pd.read_csv(infile, sep="\t", header=None, index_col=0)
    parsed = sintax.groupby(level=0).apply(parse_taxonomy, r=r, ranks=ranks, cutoff=cutoff).reset_index().drop("level_1", axis=1).set_index(0)
    parsed.index.name = "ASV"
    parsed, conf = parsed.loc[:, [x for x in parsed.columns if not x.endswith("_conf")]], parsed.loc[:, [x for x in parsed.columns if x.endswith("_conf")]]
    conf = conf.rename(columns = lambda x: x.replace("_conf", ""))
    parsed.to_csv(outfile, sep="\t")
    if conf_outfile:
        conf.to_csv(conf_outfile, sep="\t")

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("infile", type=str, help="Sintax output")
    parser.add_argument("outfile", type=str, help="Parsed outfile")
    parser.add_argument("-c", "--cutoff", type=float, default=0.8, help="Sintax cutoff (default=0.8)")
    parser.add_argument("--conf_out", type=str, help="Confidence dataframe output")
    parser.add_argument("-r", "--ranks", nargs="+", default=["kingdom", "phylum", "class", "order", "family", "genus", "species"], help="Ranks used in Sintax assignment")
    args = parser.parse_args()
    main(infile=args.infile, outfile=args.outfile, conf_outfile=args.conf_out, ranks=args.ranks, cutoff=args.cutoff)