#!/usr/bin/env python

from __future__ import unicode_literals, print_function, absolute_import
from builtins import str
import argparse as arg
import bibtexparser
import pandas as pd
import requests
import csv
import sys
import re

bare_url = "http://api.crossref.org/"

def get_bib(doi):
    """
    Parameters
    ----------

    doi: str

    Returns
    -------

    found: bool
    bib: str
    """
    url = "{}works/{}/transform/application/x-bibtex"
    url = url.format(bare_url, doi)
    r = requests.get(url)
    found = False if r.status_code != 200 else True
    bib = r.content
    bib = str(bib, "utf-8")

    return found, bib


def get_json(doi):
    """
    Parameters
    ----------
    doi: str

    Returns
    -------

    found: bool
    item: dict
        Response from crossref
    """

    url = "{}works/{}"
    url = url.format(bare_url, doi)
    r = requests.get(url)
    found = False if r.status_code != 200 else True
    item = r.json()

    return found, item


def get_bib_from_doi(doi, abbrev_journal=True, add_abstract=False):
    """
    Parameters
    ----------

    doi: str
    abbrev_journal: bool
        If True try to abbreviate the journal name
    Returns
    -------

    found: bool
    bib: str
        The bibtex string
    """
    found, bib = get_bib(doi)
    if found and abbrev_journal:

        found, item = get_json(doi)
        if found:
            abbreviated_journal = item["message"]["short-container-title"]
            if add_abstract and "abstract" in item["message"].keys():
                abstract = item["message"]["abstract"]
                bi = bibtexparser.loads(bib)
                bi.entries[0]["abstract"] = abstract
                bib = bibtexparser.dumps(bi)

            if len(abbreviated_journal) > 0:
                abbreviated_journal = abbreviated_journal[0].strip()
                bib = re.sub(
                    # r"journal = \{[^>]*?\}",
                    r"journal = \{.*\}",
                    "journal = {" + abbreviated_journal + "}",
                    bib)

    return found, bib


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(
                formatter_class=arg.ArgumentDefaultsHelpFormatter)

    #
    # Input Options
    #
    inp = parser.add_argument_group("Input Data")

    inp.add_argument('DOIs', nargs='?', default=None,
                     help='''DOIs.''')

    inp.add_argument('-i', '--input', type=str, default=None, dest='InpFile',
                     help='''Input File.''')

    #
    # Output Options
    #
    out = parser.add_argument_group("Output Options")

    out.add_argument('-o', '--output', default=None, type=str, dest='OutFile',
                     help='''Output File.''')

    args = parser.parse_args()
    Opts = vars(args)

    return Opts


if __name__ == '__main__':

    opts = options()
    if opts['InpFile'] is None:
        dois = opts['DOIs'].split()
    else:
        with open(opts['InpFile']) as f:
            dois = list(map(str.strip, f.readlines()))

    items = []
    for i, doi in enumerate(dois, start=1):
        found, item = get_bib_from_doi(doi)
        if found:
            items.append(item)
        else:
            items.append(doi)

    if opts['OutFile'] is None:
        for item in items:
            print(item)
    else:
        with open(opts['OutFile'], "w") as f:
            for item in items:
                f.write(item)
                f.write("\n%\n")
