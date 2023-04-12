import os, sys
from Bio.Seq import Seq
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                             os.pardir))
from exceptions import IterableLengthError


def ini(s, save=False, path=None, filename='rosalind_ini_1_output', ext='txt'):
    '''
    Returns tuple containig counts 'ACGT' in s 

            Parameters:
                    s (str): An input DNA string

                    save (bool): Boolean if you want to save result in a file
                    path (str): Path to either dir or file you want to save
                    filename (str): File name you want to save
                    ext (str): File extension you want to save

            Returns:
                    dna_stats (tuple): Tuple with 'ACGT' counts in s
    '''

    # assert isinstance(s, str), f'Error: type(s) = {type(s).__name__} must be str!'
    # assert len(s) < 1e3, f'Error: len(s) = {len(s)} must be < 1000!'

    if not isinstance(s, str):
        raise TypeError(f'Error: type(s) = {type(s).__name__} must be str!')
    if len(s) >= 1e3:
        raise IterableLengthError(iterable=s, sign='<', value=int(1e3))

    dna_stats = (Seq(s).count(nt) for nt in 'ACGT')

    if save:
        if path is None:
            path = os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                                f'{filename}.{ext}')
        elif os.path.isdir(path):
            path = os.path.join(path, f'{filename}.{ext}')
        with open(path, 'w') as file:
            file.write(' '.join(dna_stats))

    return dna_stats


if __name__ == '__main__':
    if len(sys.argv) == 2:
        s = sys.argv[1]
    else:
        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                               'rosalind_ini_1_dataset.txt'), 'r') as file:
            s = file.readline()

    print(s)

    print(*ini(s))
