#!/usr/env/bin python

from argparse import ArgumentParser
import pandas as pd
import polars as pl


def generate_reader(f, chunksize, nrows):
    """
    Sets up a reader with pandas. Handles both chunksize>=1 and chunksize=None

    :param f: Input file
    :param chunksize: Number of rows to read per chunk
    :param nrows: Number of total rows to read
    :return:
    """
    if nrows == 0:
        nrows = None
    r = pd.read_csv(
        f, sep="\t", index_col=0, header=0, nrows=nrows, chunksize=chunksize
    )
    if chunksize is not None:
        return r
    return [r]

