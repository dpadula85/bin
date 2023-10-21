#!/usr/bin/env python

from __future__ import unicode_literals, print_function, absolute_import
import re
import csv
import requests
import bibtexparser
import argparse as arg
from builtins import str
from impact_factor.core import Factor

from pyscopus import Scopus

bare_url = "http://api.crossref.org/"

import warnings
warnings.filterwarnings("ignore")


def options():
    '''Defines the options of the script.'''

    parser = arg.ArgumentParser(
                formatter_class=arg.ArgumentDefaultsHelpFormatter)

    #
    # Input Options
    #
    inp = parser.add_argument_group("Input Data")

    inp.add_argument(
            '-a',
            '--author',
            type=str,
            dest='Author',
            help='''Scopus author ID.''',
            default="24463809800"
        )

    inp.add_argument(
            '-k',
            '--key',
            type=str,
            required=True
            dest='APIkey',
            help='''Scopus API key.'''
        )

    #
    # Output Options
    #
    out = parser.add_argument_group("Output Options")

    out.add_argument(
            '-p',
            '--pre',
            type=str,
            dest='Prefix',
            help='''Prefix for files to save and in .bib items.''',
            default=None
        )

    out.add_argument(
            '-o',
            '--out',
            type=str,
            dest='Output',
            choices=[ "csv", "bib" ],
            help='''Format for files to save.''',
            default=[ "csv", "bib" ]
        )

    args = parser.parse_args()
    Opts = vars(args)
    if not Opts["Prefix"]:
        Opts["Prefix"] = Opts["Author"]

    return Opts


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


def get_IF(journal):
    """
    Parameters
    ----------
    journal: str

    Returns
    -------
    IF: float
    """

    engine = Factor()
    data = engine.search(journal)

    # Fix known bugs in names
    if len(data) == 0:
        journal = journal.replace(" and ", " & ")
        journal = journal.replace(" - ", "-")
        data = engine.search(journal)

    try:
        IF = data[0]['factor']
    except:
        IF = None

    return IF


if __name__ == '__main__':

    Opts = options()

    # Interact with Scopus
    scopus = Scopus(Opts["APIkey"])
    pub_df = scopus.search_author_publication(Opts["Author"])
    dois = pub_df["doi"].values.tolist()

    # Save a csv
    if "csv" in Opts["Output"]:
        df = pub_df[["scopus_id", "doi", "cover_date", "publication_name", "title"]]
        df["IF"] = df["publication_name"].apply(get_IF)
        df = df[["scopus_id", "doi", "cover_date", "IF", "publication_name", "title"]]
        df.to_csv("%s_pubs.csv" % Opts["Prefix"], index=False, quoting=csv.QUOTE_NONNUMERIC)

    # Save a bibt
    if "bib" in Opts["Output"]:
        bibitems = [ get_bib_from_doi(x)[1] for x in dois ]
        n = len(bibitems)
        bib = ""
        for i, bibitem in enumerate(bibitems):
            replacement = "@article{%s%d,\n" % (Opts["Prefix"], n - i)
            b = bibitem.split("\n")[1:]
            newbibitem = replacement + "\n".join(b)
            newbibitem += "\n%\n"
            bib += newbibitem

        with open("%s_publications.bib" % Opts["Prefix"], "w") as f:
            f.write(bib)
