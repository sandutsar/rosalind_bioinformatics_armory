import os
from Bio import SeqIO


def rcvo(file_handle, datatype='fasta', save=False, path=None, \
         filename='rosalind_rcvo_1_output', ext='txt'):
    reverse_complements = [1 if record.seq == record.seq.reverse_complement() \
                           else 0 for record in SeqIO.parse(file_handle, datatype)]

    if save:
        if path is None:
            path = os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                                f'{filename}.{ext}')
        elif os.path.isdir(path):
            path = os.path.join(path, f'{filename}.{ext}')
        with open(path, 'w') as file:
            file.write(str(sum(reverse_complements)))

    return sum(reverse_complements)


if __name__ == '__main__':
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                            'rosalind_rcvo_1_dataset.txt'), 'r') as file:
        print(rcvo(file, save=True))
