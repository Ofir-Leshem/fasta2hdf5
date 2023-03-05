# fasta2hdf5

### Description
<p>
This project converts FASTA files into HDF5 files containing mean currents values.<br>
The FASTA files are divided into k-mers (k=6), and matched with a mean current estimate.
</p>

### Quickstart

#### Update config file
<p>
update SRC_DIR to hold the path in which the FASTA files are saved.<br>
update DST_DIR to hold the path to which the HDF5 files will be saved.<br>
</p>

#### Script
<p>
The fasta2hdf5.py script does the translation.<br>
The 6mers.txt file should be saved in the same directory as the script.
</p>

### Credits
<p>
Our work is based on the following:<br>
<br>
k-mer models were taken from: <br>
https://github.com/langjidong/ONT-Nanopore-kmer_models
<br>
<br>
SquiggleFilter:<br>
https://github.com/TimD1/SquiggleFilter<br>
Dunn, Tim, Harisankar Sadasivan, Jack Wadden, Kush Goliya, Kuan-Yu Chen, David Blaauw, Reetuparna Das, and Satish Narayanasamy. "SquiggleFilter: An Accelerator for Portable Virus Detection." In MICRO-54: 54th Annual IEEE/ACM International Symposium on Microarchitecture, pp. 535-549. 2021.
<br>
</p>


### Notes
<p>
The project currently ignores 'N' bases altogether, but support will be added.
</p>

