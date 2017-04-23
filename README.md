# transporter_annotation
#### Annotation of membrane bound transporter proteins from BLASTP, TIGRfam, Pfam and TMHMM

Written by John P. McCrow (J. Craig Venter Institute, La Jolla, CA 92037).
Adapted from previous work by Johnathan Badger (https://github.com/jhbadger/scripts)
The database was originally developed by Qinghu Ren and Ian Paulsen, and is based on TransportDB (Qinghu Ren, et al., TransportDB: a comprehensive database resource for cytoplasmic membrane transport systems and outer membrane channels. Nucleic Acids Research 35: 274-279, 2007)

### *Under active development. Check back later for a stable version...*

Usage
-----

| File | Description |
|------|-------------|
| find_transporters.py | Main executable |
| | |
| trans_db.sqlite3 | Database of transporter annotations |
| trans_models.hmm | Pfam and TIGRfam HMMs |
| trans_peps.fa | Peptides FASTA |


```
find_transporters v0.1 (Apr 13, 2016)
Membrane-bound transporter prediction

Usage: ./find_transporters.py (options)
   -f file        : FASTA file (required if missing any of -b,-p,-t)
   -b file        : BLASTP output file (required if missing -f)
   -p file        : HMM output file (required if missing -f)
   -t file        : TMHMM output file (required if missing -f)

   -c, --cpus int : parallel processes for BLASTP, HMMER, TMHMM (default: 3)
   -h, --help     : help
   -l             : extended output format
   -m num         : minimum evidence score to output (default: 2)
   -n             : only accept input files (-b,-p,-t); do NOT run processes
   -o file        : output file (default: stdout)
   -s             : brief output format (default)
   -v, --verbose  : more information to stderr

FASTA file or output from all BLAST, HMM, TMHMM required.
```

Installation
------------

Install all dependencies and then run the make utility to index sequence and hmm files:
```
make
```

Output Format
-------------

Output is a tab-delimited table with the following columns:
peptide id
transporter family
transporter substrate
annotation score

If extended output format is selected using the -l parameter, then the output fields are:
peptide id
transporter family
transporter sub-family
transporter substrate
annotation score
transporter classification (tc) id
trans-membrane helices
trans-membrane length
trans-membrane topology

Dependencies
------------

* Python (https://www.python.org/downloads/)
* Sqlite3 (https://sqlite.org/download.html)
* Hmmer3 (http://hmmer.org/download.html)
* NCBI-Blast (https://blast.ncbi.nlm.nih.gov/Blast.cgi)
* TMHMM (http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm)

*The last 3 dependencies are not required to run the program, as long as output from all 3 are present.*
