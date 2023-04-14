import os, sys
from Bio import Entrez
Entrez.email = 'sandutsar@yahoo.com'


def gbk(organism, mindate, maxdate, db='nucleotide', datetype='pdat', \
        save=False, path=None, filename='rosalind_gbk_1_output', ext='txt'):
    records_count = Entrez.read(Entrez.esearch(db=db, term=f'{organism}[Organism]', \
                                               datetype=datetype, mindate=mindate, \
                                                maxdate=maxdate))['Count']

    if save:
        if path is None:
            path = os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                                f'{filename}.{ext}')
        elif os.path.isdir(path):
            path = os.path.join(path, f'{filename}.{ext}')
        with open(path, 'w') as file:
            file.write(str(records_count))

    return records_count


if __name__ == '__main__':
    if len(sys.argv) == 4:
        organism, mindate, maxdate = sys.argv[1:]
    else:
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                               'rosalind_gbk_1_dataset.txt'), 'r') as file:
            organism, mindate, maxdate = file.readlines()
   
    print(organism, mindate, maxdate)

    print(gbk(organism, mindate, maxdate, save=True))
