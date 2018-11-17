
from bs4 import BeautifulSoup
import requests
import argparse
import os
import socket

class PeptideCutterLoop():
    """
    this class takes a table from website, containing following columns:
    """
    def __init__(self, table_pep_id_coords):
        self.table_pep_id_coords = table_pep_id_coords
        self.outFile = open(self.table_pep_id_coords + "_parsed.pepcutter", "w")
        self.run()


    def getResponse(self,id):
        """
        Send post request to the server
        :param id: protein id
        :return: response object
        """
        headers = {'User-Agent': 'Mozilla/5.0'}
        payload = {'protein': id, 'enzyme_number': 'all_enzymes', 'cleave_number': 'all', 'alphtable': 'alphtable'}

        session = requests.Session()
        response = session.post("https://web.expasy.org/cgi-bin/peptide_cutter/peptidecutter.pl/", headers=headers,
                                data=payload)
        return response

    def parseResponse(self,response, id):
        """
        parsing response from the server
        :param response: response object from self.getResponse
        :param id: proteins_id
        :return:
        write parsing results in separate file "./PC_tmp_files/table_{}.txt".format(id)

        """
        soup = BeautifulSoup(response.content, "lxml")
        out = open("./PC_tmp_files/table_{}.txt".format(id), "w")
        table = soup.findAll("table", class_="proteomics2")[0]
        # print(table)
        cnt = 0
        for real_rows in table.findAll("tr"):
            cnt += 1
            # print("<<<<<<<<<", cnt)
            # print(real_rows)
            td = [i.text.replace("\n", "") for i in real_rows.findAll("td")]
            out.write("\t".join(td) + "\n")
        out.close()
        # time.sleep(5)

    def getOverlap(self, pep_coords,cleavege_file):
        """
        compare coordinates of peptide and cleaveg sites
        :param pep_coords:
        :param cleavege_file: file from
        :return:
        """
        pep_coords = [coords.split("-") for coords in pep_coords.split(", ")]
        int_coords = []
        for coords in pep_coords:
            print(coords)
            int_coords.append(sorted([int(coords[0]), int(coords[1])]))

        with open(cleavege_file) as infile:
            enzymes_start = []
            enzymes_stop = []
            for lines in infile:
                sp = lines.rstrip().split("\t")
                if len(sp)>1:
                    enzyme = sp[0]
                    sites = [int(sites) for sites in sp[2].split(" ")]

                    for site in sites:
                        for pep_c in int_coords:
                            if site == pep_c[0]:
                                enzymes_start.append(enzyme)
                            elif site == pep_c[1]:
                                enzymes_stop.append(enzyme)
            if enzymes_start:
                enzymes_start = ",".join([en for en in set(enzymes_start)])
            else:
                enzymes_start = "-"

            if enzymes_stop:
                enzymes_stop = ",".join([en for en in set(enzymes_stop)])
            else:
                enzymes_stop = "-"

            self.outFile.write(enzymes_start +"\t" + enzymes_stop + "\n")


    def run(self):
        with open(self.table_pep_id_coords) as infile:
            self.outFile.write("\t".join(['protein_id', 'peptide_seq' ,"Name of enzyme overlap START", "Name of enzyme overlap END"]) + "\n")
            for i, lines in enumerate(infile):
                sp = lines.rstrip().split("\t")
                id = sp[0]
                print(i, id)
                peptide_seq = sp[2]
                peptides_coords = sp[1]
                resp = self.getResponse(id)
                self.parseResponse(resp, id)
                self.outFile.write(id + "\t" + peptide_seq + "\t")
                self.getOverlap(peptides_coords,"./PC_tmp_files/table_{}.txt".format(id))

        self.outFile.close()

if __name__ == '__main__':
    def is_connected(hostname):
        try:
            # see if we can resolve the host name -- tells us if there is
            # a DNS listening
            host = socket.gethostbyname(hostname)
            # connect to the host -- tells us if the host is actually
            # reachable
            s = socket.create_connection((host, 80), 2)
            return True
        except:
            pass
        return False

    parser = argparse.ArgumentParser(description='This script is useful for overlapping the protease sites cleavage and sites of MS peptide origins.')
    parser.add_argument("table", help='a tab separated table with first two columns are (!no header!): 	\n'
                                      '1. Protein ID (!recognizable by PeptideCutter) \n \
	                                   2. Sites for MS peptide origin (e.g. 23-35) \n \
                                       3. Peptide sequence')


    args = parser.parse_args()
    if not os.path.isfile(args.table):
        print("No files were found! {}".format(args.table))
    else:
        print("File {} found...".format(args.table))

    if not is_connected("https://web.expasy.org/peptide_cutter/"):
        print("PeptideCutter does not respond")
    else:
        "Peptide cutter server is ON..."

    if not os.path.isdir("PC_tmp_files"):
        os.makedirs("PC_tmp_files")
    print("*****************RUN*****************")
    PeptideCutterLoop(args.table)



