# from pyteomics import parser
# products = parser.cleave('AKAKBK', parser.expasy_rules['trypsin'], 0)
# print(products)
from Bio import SeqUtils
from bs4 import BeautifulSoup
import requests
from collections import defaultdict
import re

def makeRE(matrix_dictionary):
    re_pattern = ''
    for positios in range(8):
        list_aa = ""
        for letters in matrix_dictionary:
            fq_aa_in_position = matrix_dictionary[letters][positios]
            if fq_aa_in_position == 0:
                list_aa += letters
        if list_aa == "":
            re_pattern += "."
        else:
            re_pattern += "[^{}]".format(list_aa)
    print(re_pattern)
    return re_pattern


output_file = open("RE_file.tab",'w')
def parseEnzyme_html(html, enzyme):
    result = requests.get(html).content
    soup = BeautifulSoup(result,'lxml')
    tab = soup.find("table", {'summary':'matrix'})
    matrix_dictionary = {}
    if tab:
        for i,rows in enumerate(tab.findAll('tr')):
            values = [row_val.text for row_val in rows.findAll("td")]
            if values:
                aa_leter = SeqUtils.IUPACData.protein_letters_3to1[values[0]]
                frequencies = [int(i) for i in values[1:]]
                matrix_dictionary[aa_leter] = frequencies
        print(matrix_dictionary)
        re_patters = makeRE(matrix_dictionary)
        output_file.write(enzyme + "\t" + re_patters + "\n")
        return True
    else:
        return False



def parseMEROPS():
    cnt_no_matrix = 0
    cnt_matrix = 0
    root = "https://www.ebi.ac.uk"
    html = "https://www.ebi.ac.uk/merops/cgi-bin/speccards?sp=sp002565;type=peptidase"
    result = requests.get(html).content
    soup = BeautifulSoup(result, 'lxml')
    tab = soup.find("table", {'id': 'details'})
    tab_rows = tab.findAll('tr')
    for i, rows in enumerate(tab_rows):
        print(rows)
        cols = rows.findAll('td')
        if cols:
            enzyme = cols[3].text
            href = root + cols[2].find('a')['href']
            if not parseEnzyme_html(href, enzyme):
                cnt_no_matrix += 1
            else:
                cnt_matrix += 1

    print("Number of enzymes with NO cleavege matrix: ", cnt_no_matrix)
    print("Number of enzymes with cleavege matrix: ", cnt_matrix)

# html = "https://www.ebi.ac.uk/merops/cgi-bin/pepsum?id=A01.020"
# parseEnzyme_html(html)
parseMEROPS()

output_file.close()
