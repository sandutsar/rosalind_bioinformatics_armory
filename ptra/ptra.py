import os, sys
from Bio.Seq import translate
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                             os.pardir))
from exceptions import IterableLengthError


def ptra(dna, protein, save=False, path=None, \
         filename='rosalind_ptra_1_output', ext='txt'):
    # assert len(dna) <= 1e4, f'Error: len(dna) = {len(dna)} must be <= 10000'
    if len(dna) > 1e4:
        raise IterableLengthError(iterable=dna, sign='<=', value=int(1e4))
    
    tables = []
    for table in [x for x in range(1, 16)]:
        try:
            if translate(dna, table, to_stop=True) == protein:
                tables += [table]
        except KeyError:
            continue

    if save:
        if path is None:
            path = os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                                f'{filename}.{ext}')
        elif os.path.isdir(path):
            path = os.path.join(path, f'{filename}.{ext}')
        with open(path, 'w') as file:
            file.write(str(tables[0]))
    
    return tables[0]


if __name__ == '__main__':
    if len(sys.argv) == 3:
        dna = sys.argv[1]
        protein = sys.argv[2]
    else:
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                            'rosalind_ptra_1_dataset.txt'), 'r') as file:
            dna = file.readline()[:-1]
            protein = file.readline()[:-1]

    print(dna, protein)

    print(ptra(dna, protein, save=True))
