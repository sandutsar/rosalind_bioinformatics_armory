import os, sys
from Bio import Entrez, SeqIO
Entrez.email = 'sandutsar@yahoo.com'
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                             os.pardir))
from exceptions import IterableLengthError


def frmt(ids, db='nucleotide', rettype='fasta', save=False, path=None, \
         filename='rosalind_frmt_1_output', ext='txt'):
    if len(ids) > 10:
        raise IterableLengthError(iterable=ids, sign='<=', value=10)
    
    records = list(SeqIO.parse(Entrez.efetch(db=db, id=ids, \
                                              rettype=rettype), rettype))
    shortest_string_record = records[[len(record.seq) for record in \
        records].index(min([len(record.seq) for record in records]))]
    shortest_string = [shortest_string_record.description, \
                       shortest_string_record.seq]
    
    if save:
        if path is None:
            path = os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                                f'{filename}.{ext}')
        elif os.path.isdir(path):
            path = os.path.join(path, f'{filename}.{ext}')
        with open(path, 'w') as file:
            SeqIO.write(shortest_string_record, file, rettype)
    
    return shortest_string


if __name__ == '__main__':
    if len(sys.argv) > 1 and len(sys.argv) < 12:
        ids = sys.argv[1:]
    else:
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                                'rosalind_frmt_1_dataset.txt'), 'r') as file:
            ids = file.readline().split()

    print(*ids)

    print(*frmt(ids, save=True))
