import math
import os
import random

import h5py
import numpy as np
from Bio import SeqIO
import configparser

from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove

READ_SIZE = 1000000
NUMBER_OF_READS = 1


def rand(sequences):
    bases = ['A', 'C', 'G', 'T']
    rand_sequences = []
    for seq in sequences:
        for i in range(len(seq)):
            if seq[i] != 'A' and seq[i] != 'C' and seq[i] != 'G' and seq[i] != 'T':
                seq = seq[:i] + bases[random.randint(0, 3)] + seq[i + 1:]
        rand_sequences.append(seq)
    return rand_sequences


def readFastaSequence(fasta_filepath):
    fasta_sequences = SeqIO.parse(open(fasta_filepath), 'fasta')
    name, sequence = [], []
    for fasta in fasta_sequences:
        # sequence.append(str(fasta.seq))
        rem = len(fasta) % READ_SIZE
        if rem == 0:
            split = np.array_split(fasta.seq, len(fasta.seq) / READ_SIZE)
        elif len(fasta) >= READ_SIZE:
            split = np.array_split(fasta.seq[:-rem], len(fasta.seq) / READ_SIZE)
        else:
            split = []

        for s in split:
            if len(s) >= READ_SIZE:
                name.append(fasta.id)
                sequence.append(''.join(s))

    return rand(sequence)


def build_6mers_dict():
    d = {}
    with open("6mers.txt", 'r') as f:
        for line in f:
            val = line.split()
            d[val[0]] = float(val[1])
    return d


def translate(sequences, kmers_dict, kmer_size=6):
    normalized = []
    for seq in sequences:
        l = np.zeros(len(seq))
        cnt = 0
        for k in range(len(seq) - kmer_size):
            if 'N' in seq[k:k + kmer_size]:
                cnt += 1
            else:
                l[k - cnt] = kmers_dict[seq[k:k + kmer_size]]
        for i in range(len(l)):
            if math.isnan(l[i]) or l[i] == 0:
                np.delete(l, i)
        normalized.append(discrete_normalize(l * 100))
    return normalized


# Code by https://github.com/TimD1/SquiggleFilter/blob/master/sdtw_analysis.ipynb
def discrete_normalize(seq, bits=8, minval=-4, maxval=4):
    """ Approximate normalization which converts signal to integer of desired precision. """
    mean = int(np.mean(seq))
    mean_avg_dev = int(np.mean(np.abs(seq - mean)))
    norm_seq = (seq - mean) / mean_avg_dev
    norm_seq[norm_seq < minval] = minval  # threshold
    norm_seq[norm_seq > maxval] = maxval
    norm_seq = ((norm_seq - minval) * (2 ** bits / (maxval - minval))).astype(int)
    return norm_seq


def create_hdf5(values, hdf5_filepathpath):
    if len(values) >= NUMBER_OF_READS:
        n_samples = len(values)
        rand_gen = np.random.RandomState(0)
        indices = np.arange(n_samples)
        # rand_gen.shuffle(indices)
        indices = indices[:NUMBER_OF_READS]
        values = [values[i] for i in indices]
        f = h5py.File(hdf5_filepathpath, 'w')
        for i in range(len(values)):
            s = "read_" + str(i)
            f.create_dataset(s, data=values[i])


def fasta2hdf5(fasta_path, hdf5_path):
    kmers_dict = build_6mers_dict()
    for fn in os.listdir(fasta_path):
        if os.path.isfile(os.path.join(fasta_path, fn)) and (
                os.path.splitext(fn)[-1].lower() == ".fasta" or os.path.splitext(fn)[-1].lower() == ".fna"):
            seq = readFastaSequence(os.path.join(fasta_path, fn))
            currents = translate(seq, kmers_dict)
            create_hdf5(currents, os.path.join(hdf5_path, os.path.splitext(fn)[0] + ".hdf5"))


if __name__ == '__main__':
    # replace(r"C:\Users\ofirl\PycharmProjects\pythonProject\pythonProject\NCBI_data\e_coli\GCF_946909605.1_NILS35_genomic.fna")
    config = configparser.ConfigParser()
    config.read("config.ini")
    src_path = config.get('INFO', 'SRC_DIR')
    dst_path = config.get('INFO', 'DST_DIR')
    fasta2hdf5(src_path, dst_path)
