import urllib.request
import json
import argparse
import re

def msgfToUnipept(peptides_in):

    """
    Takes the list of peptides and, cleaves them into Unipept format (0 missed cleavages,
    cleaves after K or R except followed by P
    :param peptides_in: text file of peptides, one peptide per line
    :return: list of peptides in Unipept format
    """

    peptides = list()
    trypsin = lambda peptide: re.sub(r'(?<=[RK])(?=[^P])', '\n', peptide, re.DOTALL).split()
    with open(peptides_in, "r") as pep_in:
        for peptide in pep_in:
            peptides += trypsin(peptide)

    return peptides



def generateGetRequest(peptides):

    """
    Given a list of Unipept peptides, it generates a list of Unipept Get Requests
    :param peptides: Unipept peptides
    :return: list of get requests for Unipept
    """

    request_list = list()

    char_count = 74  # number of chars minimum per request
    request = "http://api.unipept.ugent.be/api/v1/peptinfo.json?"
    for peptide in peptides:
        if (len(peptide) + 9 + char_count) < 2048:
            request += "input[]={}&".format(peptide)
            char_count += 9 + len(peptide)
        else:
            request += "equate_il=true&extra=true&names=true"
            request_list.append(request)
            request = "http://api.unipept.ugent.be/api/v1/peptinfo.json?" + "input[]={}&".format(peptide)
            char_count = 74 + 9 + len(peptide)
    if len(request) != 49:  # first part of request
        request += "equate_il=true&extra=true"
        request_list.append(request)

    return request_list


def getInfoFromUnipept(request_list, result_file):

    """
    Send all requests, get for each peptide the phylum, family, genus and collection of EC-numbers
    :param request_list: list of Get Requests
    :param result_file: csv file with Unipept info (phylum, family, genus and collection of EC-numbers)
    :return: None
    """

    unipept_list = list()
    for element in request_list:
        unipept_json = urllib.request.urlopen(element).read()
        unipept_list.append(json.loads(unipept_json.decode('utf-8')))

    tmp = open(result_file, "w")
    tmp.close()

    with open(result_file, "a") as f_out:
        for response in unipept_list:
            for element in response:
                try:
                    kingdom_name = element["kingdom_name"].strip()
                except KeyError:
                    kingdom_name = ""
                try:
                    phylum_name = element["phylum_name"].strip()
                except KeyError:
                    phylum_name = ""
                try:
                    class_name = element["class_name"].strip()
                except KeyError:
                    class_name = ""
                try:
                    order_name = element["order_name"].strip()
                except KeyError:
                    order_name = ""
                try:
                    family_name = element["family_name"].strip()
                except KeyError:
                    family_name = ""
                try:
                    genus_name = element["genus_name"].strip()
                except KeyError:
                    genus_name = ""
                try:
                    species_name = element["species_name"].strip()
                except KeyError:
                    species_name = ""
                try:
                    ec_list = list()
                    for i in element['ec']:
                        ec_list.append(i['ec_number'].strip())
                    ec = ";".join(ec_list)
                except IndexError:
                    ec = ""

                print("{},{},{},{},{},{},{},{},{}".format(element["peptide"].strip(), kingdom_name, phylum_name, class_name, order_name, family_name, genus_name, species_name, ec), file=f_out)


if __name__ == "__main__":

    print("""
Taxonomical and Functional analysis using Unipept API
Equate I/L = True; Filter duplicate peptides = False
    """)

    # parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="peptides in text format, one peptide per line", type=str)
    parser.add_argument("-o", help="output file in csv format", type=str)
    args = parser.parse_args()

    # from peptide list from MS-GF+ to Unipept peptides
    peptides = msgfToUnipept(args.i)

    # Generate list of Get Requests for Unipept
    request_list = generateGetRequest(peptides)

    # Send all requests, get for each peptide the phylum, family, genus and collection of EC-numbers
    getInfoFromUnipept(request_list, args.o)



