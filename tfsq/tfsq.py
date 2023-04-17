import os, sys
from Bio import SeqIO


def tfsq(file_handle, input_type='fastq', output_type='fasta', save=True, path=None, \
         filename='rosalind_tfsq_1_output', ext='txt'):
    data = SeqIO.parse(file_handle, input_type)

    if save:
        if path is None:
            path = os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                                f'{filename}.{ext}')
        elif os.path.isdir(path):
            path = os.path.join(path, f'{filename}.{ext}')
        with open(path, 'w') as file:
            SeqIO.write(data, file, output_type)


if __name__ == '__main__':
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                                'rosalind_tfsq_1_dataset.txt'), 'r') as file:
        tfsq(file)
