# BLAST Discrimination
Using BLAST on a FASTA file to determine taxonomy discrimination.

# Prerequresits
Requires the NCBI BLAST+ tools to be installed and on your execution path. To install on Ubuntu:

```bash
sudo apt-get install ncbi-blast+
```

The following python packages are required to run the script:

```bash
python3 -m pip install pandas
python3 -m pip install biopython
```

# How to run
```bash
python3 discrimination.py --fasta <fasta-file> --meta <excel-meta-data>
```

## Arguments
 --fasta
 
           Fasta file that has identifiers formatted as below

           ><unique-id> <genus> <species>
           <..sequence..>
           ><unique-id> <genus> <species>
           <..sequence..>

           The full fasta file is used to build the BLAST database. Only
           sequences for species that have >= 2 sequences will be used within
           the discrimination analysis.

 --meta
 
           Excel meta-data that has the following column header names
           (at least). Note that "order" is optional but the rest are
           mandatory.

           Accepted_Name, [order], family, genus, species

           The Accepted_Name is used to match against the <genus> <species>
           from the fasta file (see above).

## Output
* A blast database made from the input fasta file "blastdb.*"
* A FASTA file formatted to build the blast database "blast_discrimination_fasta.fas"
* Working file "output.csv" that is a single output from a BLAST query
* Working file "query.fas" that is a single input to the BLAST query
* Discrimination results "blast_discrimination_fasta.fas"
