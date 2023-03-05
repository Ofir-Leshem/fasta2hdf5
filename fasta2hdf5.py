import math
import os

import h5py
import numpy as np
from Bio import SeqIO


def readFastaSequence(fasta_filepath):
    fasta_sequences = SeqIO.parse(open(fasta_filepath), 'fasta')
    name, sequence = 0, 0
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
    return sequence


def build_6mers_dict():
    d = {}
    with open("6mers.txt", 'r') as f:
        for line in f:
            val = line.split()
            d[val[0]] = float(val[1])
    return d


def translate(seq, kmers_dict, kmer_size=6):
    l = np.zeros(len(seq))
    cnt = 0
    for k in range(len(seq) - kmer_size):
        if 'N' in seq[k:k + kmer_size]:
            cnt += 1
        elif 'N' in seq[k:k + kmer_size]:
            l[k-cnt] = find_avg(seq[k:k + kmer_size], kmers_dict)
        else:
            l[k-cnt] = kmers_dict[seq[k:k + kmer_size]]
    for i in range(len(l)):
        if math.isnan(l[i]):
            np.delete(l, i)
    return discrete_normalize(l*100)


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


def find_avg(seq, kmers_dict):
    bases = ['A', 'C', 'G', 'T']
    base_seq = ['A', 'A', 'A', 'A', 'A', 'A']
    loc = []
    curr_base_v = []
    for i in range(len(seq)):
        base_seq[i] = seq[i]
        if seq[i] == 'N':
            loc.append(i)
            base_seq[i] = bases[0]
            curr_base_v.append(0)
    last = base_seq
    curr = last
    k = len(loc)-1
    current_sum = kmers_dict[''.join(curr)]
    for i in range(np.power(len(bases), len(loc))-1):
        k = len(loc)-1
        while k != -1:
            if curr_base_v[k] < len(bases)-1:
                curr_base_v[k] += 1
                curr[loc[k]] = bases[curr_base_v[k]]
                k = -1
            elif k != 0 and curr_base_v[k] == len(bases)-1:
                curr_base_v[k] = 0
                curr[loc[k]] = bases[0]
                k -= 1
            else:
                return np.divide(current_sum,np.power(len(bases), len(loc)))

        current_sum += kmers_dict[''.join(curr)]
        last = curr


def create_hdf5(values, hdf5_filepathpath):
    f = h5py.File(hdf5_filepathpath, 'w')
    dset = f.create_dataset("dataset", data=values)


def fasta2hdf5(fasta_path, hdf5_path):
    kmers_dict = build_6mers_dict()
    for fn in os.listdir(fasta_path):
        if os.path.isfile(os.path.join(fasta_path, fn)) and os.path.splitext(fn)[-1].lower() == ".fasta":
            seq = readFastaSequence(os.path.join(fasta_path, fn))
            currents = translate(seq, kmers_dict)
            create_hdf5(currents, os.path.join(hdf5_path, os.path.splitext(fn)[0]+".hdf5"))


if __name__ == '__main__':
    src_path = ""
    dst_path = ""
    fasta2hdf5(src_path, dst_path)
