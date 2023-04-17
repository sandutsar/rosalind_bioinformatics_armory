import os
from statistics import mean
from Bio import SeqIO


def phre(file_handle, threshold, datatype='fastq', save=False, path=None, \
         filename='rosalind_phre_1_output', ext='txt'):
    avg_qual_below_thres = sum([1 if mean(record.letter_annotations['phred_quality']) \
                                < threshold else 0 for record in SeqIO.parse(file_handle, datatype)])

    if save:
        if path is None:
            path = os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                                f'{filename}.{ext}')
        elif os.path.isdir(path):
            path = os.path.join(path, f'{filename}.{ext}')
        with open(path, 'w') as file:
            file.write(str(avg_qual_below_thres))

    return avg_qual_below_thres


if __name__ == '__main__':
    with open(os.path.join(os.path.dirname(os.path.realpath(__file__)), \
                                'rosalind_phre_1_dataset.txt'), 'r') as file:
        threshold = int(file.readline())
        print(phre(file, threshold, save=True))
