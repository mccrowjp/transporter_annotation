#!/usr/bin/env python
#
# find_transporters v0.1 (Apr 13, 2016)
# Membrane-bound transporter prediction
#

import sys, re, os, getopt
import bz2, gzip, time
import multiprocessing, subprocess, signal
import sqlite3

prog_path = os.path.realpath(sys.argv[0])
prog_dir = os.path.dirname(prog_path)

trans_db_file = os.path.join(prog_dir, 'trans_db.sqlite3')
trans_hmm_file = os.path.join(prog_dir, 'trans_models.hmm')
trans_blast_file = os.path.join(prog_dir, 'trans_peps.fa')

# transporter score weighting
weight_cog_positive = 1
weight_cog_negative = -2
weight_tc = 1
weight_pfam = 1
weight_tm_helices = 1

# thresholds
min_tm_helices = 3
max_blast_evalue = 1e-5
max_hmm_evalue = 1e-4

# global defaults
cpus = 3
min_orf_score = 2.0
verbose = False

run_type_queue = []

class Run_type:
    blast = 1
    hmm = 2
    tmhmm = 3

class Hmm_format:
    unknown = 0
    tblout = 1
    htab = 2

class Run_list_element:
    def __init__(self, type, infile, outfile):
        self.type = type
        self.infile = infile
        self.outfile = outfile

class Transporter:
    def __init__(self, family, substrate):
        self.family = family
        self.substrate = substrate

class Hmm_info:
    def __init__(self, tc, family, subfamily):
        self.tc = tc
        self.family = family
        self.subfamily = subfamily

class Cog_info(Hmm_info):
    def __init__(self, tc, family, subfamily, substrate):
        Hmm_info.__init__(self, tc, family, subfamily)
        self.substrate = substrate

class Orf_hits:
    id = ""
    pfam_id = ""
    pfam_evalue = 0
    tm_length = 0
    tm_helices = 0
    tm_topo = ""
    blast_tc = ""
    blast_evalue = 0
    blast_cog = ""
    blast_cog_evalue = 0
    score = 0
    
    def __init__(self, id):
        self.id = id

dict_transporter = {}
dict_hmm_info = {}
dict_cog_info = {}
dict_tc_seguid = {}
dict_cog_seguid = {}
dict_fp_cog = {}
dict_orf = {}

def xprint(s):
        sys.stderr.write(str(s) + '\n')
        
def run_all_parallel():
    if cpus > 1:
        pool = multiprocessing.Pool(cpus)
        pool.map(run_each, run_type_queue)
        try:
            pool.get(1)
        except KeyboardInterrupt:
            pool.terminate()
            pool.join()
        else:
            pool.close()
            pool.join()

    else:
        for rlt in run_type_queue:
            run_each(rlt)

def run_each(rlt):
    try:
        if rlt.type == Run_type.blast:
            run_blast(rlt.infile, rlt.outfile)
        elif rlt.type == Run_type.hmm:
            run_hmm(rlt.infile, rlt.outfile)
        elif rlt.type == Run_type.tmhmm:
            run_tmhmm(rlt.infile, rlt.outfile)
    except KeyboardInterrupt:
        pass


def run_blast(infile, outfile):
    if verbose:
        xprint("Running BLAST: " + infile + " -> " + outfile)
    
    rc = subprocess.call("blastall -p blastp -d " + trans_blast_file + " -i " + infile + " -o " + outfile + " -m8 -e " + str(max_blast_evalue) + " 2>/dev/null", shell=True)
    xprint("[run_blast]: blastall = " + str(rc))

def run_hmm(infile, outfile):
    if verbose:
        xprint("Running HMM: " + infile + " -> " + outfile)

    rc = subprocess.call("hmmscan --tblout " + outfile + " " + trans_hmm_file + " " + infile + " &>/dev/null", shell=True)
    xprint("[run_hmm]: hmmscan = " + str(rc))

def run_tmhmm(infile, outfile):
    if verbose:
        xprint("Running TMHMM: " + infile + " -> " + outfile)

    rc = subprocess.call("tmhmm " + infile + " > " + outfile + " 2>/dev/null", shell=True)
    xprint("[run_tmhmm]: tmhmm = " + str(rc))

def read_trans_db():
    global dict_transporter
    global dict_hmm_info
    global dict_cog_info
    global dict_tc_seguid
    global dict_cog_seguid
    global dict_fp_cog
    
    if verbose:
        xprint("Reading DB: " + trans_db_file)

    transdb_con = sqlite3.connect(trans_db_file)
    transdb_cursor = transdb_con.cursor()

    for row in transdb_cursor.execute('SELECT tc, family, substrate FROM transporter'):
        dict_transporter[row[0]] = Transporter(row[1], row[2])

    for row in transdb_cursor.execute('SELECT name, tc, family, subfamily FROM hmm_info'):
        dict_hmm_info[row[0]] = Hmm_info(row[1], row[2], row[3])

    for row in transdb_cursor.execute('SELECT name, tc, family, subfamily, substrate FROM cog_info'):
        dict_cog_info[row[0]] = Cog_info(row[1], row[2], row[3], row[4])

    for row in transdb_cursor.execute('SELECT seguid, tc FROM tc_seguid'):
        dict_tc_seguid[row[0]] = row[1]

    for row in transdb_cursor.execute('SELECT seguid, cog FROM cog_seguid'):
        dict_cog_seguid[row[0]] = row[1]

    for row in transdb_cursor.execute('SELECT name FROM fp_cog'):
        dict_fp_cog[row[0]] = 1


def read_blast(blast_file):
    global dict_orf
    in_handle = happyfile.hopen_or_else(blast_file)

    if verbose:
        xprint("Reading BLAST file: " + blast_file)

    while 1:
        line = in_handle.readline()
        if not line:
            break
        line = line.rstrip()
        qid, sid, pid, len, mm, go, qs, qe, ss, se, estr, bs = line.split("\t")[:12]

        evalue = float(estr)
    
        if evalue <= max_blast_evalue:
            if sid in dict_tc_seguid:
                orf = dict_orf.get(qid, Orf_hits(qid))
                if (not orf.blast_tc) or evalue < orf.blast_evalue:
                    orf.blast_tc = dict_tc_seguid[sid]
                    orf.blast_evalue = evalue
                    dict_orf[orf.id] = orf
            if sid in dict_cog_seguid:
                orf = dict_orf.get(qid, Orf_hits(qid))
                if (not orf.blast_cog) or evalue < orf.blast_cog_evalue:
                    orf.blast_cog = dict_cog_seguid[sid]
                    orf.blast_cog_evalue = evalue
                    dict_orf[orf.id] = orf

def read_hmms(hmm_file):
    global dict_orf
    in_handle = happyfile.hopen_or_else(hmm_file)
    
    if verbose:
        xprint("Reading HMM file: " + hmm_file)

    hmm_format = Hmm_format.unknown
    line_num = 0
    while 1:
        line = in_handle.readline()
        line_num += 1
        if not line:
            break
        line = line.rstrip()

        if line_num == 1:
            if re.match('\#\s+\-+ full', line):
                hmm_format = Hmm_format.tblout
            elif len(line.split("\t")) == 21:
                hmm_format = Hmm_format.htab

        if not re.match('\#', line):
            hid = ""
            qid = ""
            estr = 'Inf'
            if hmm_format == Hmm_format.tblout:
                cols = re.split('\s+', line)
                hid, hidfull, qid, acc, estr = cols[:5]
            elif hmm_format == Hmm_format.htab:
                cols = re.split('\t', line)
                hid = cols[0]
                qid = cols[5]
                estr = cols[19]
                        
            hid = re.sub('\.\d+$', '', hid)
            evalue = float(estr)
            if evalue <= max_hmm_evalue:
                if hid in dict_hmm_info:
                    orf = dict_orf.get(qid, Orf_hits(qid))
                    if (not orf.pfam_id) or evalue < orf.pfam_evalue:
                        orf.pfam_id = hid
                        orf.pfam_evalue = evalue
                        dict_orf[orf.id] = orf


def read_tmhmm(tmhmm_file):
    global dict_orf
    in_handle = happyfile.hopen_or_else(tmhmm_file)
    
    if verbose:
        xprint("Reading TMHMM file: " + tmhmm_file)
    
    while 1:
        line = in_handle.readline()
        if not line:
            break
        line = line.rstrip()

        m1 = re.match('^\# (\S+) Length: (\d+)', line)
        if m1:
            orf = dict_orf.get(m1.group(1), Orf_hits(m1.group(1)))
            orf.tm_length = int(m1.group(2))
            dict_orf[orf.id] = orf

        m2 = re.match('^\# (\S+) Number of predicted TMHs:\s+(\d+)', line)
        if m2:
            orf = dict_orf.get(m2.group(1), Orf_hits(m2.group(1)))
            orf.tm_helices = int(m2.group(2))
            dict_orf[orf.id] = orf

        m3 = re.match('^(\S+)\t\S+\t(inside|outside|TMhelix)\t(.+)$', line)
        if m3:
            orf = dict_orf.get(m3.group(1), Orf_hits(m3.group(1)))
            if m3.group(2) == 'inside':
                orf.tm_topo += "i"
            elif m3.group(2) == 'outside':
                orf.tm_topo += "o"
            elif m3.group(2) == 'TMhelix':
                m4 = re.search('\s*(\d+)\s+(\d+)', m3.group(3))
                orf.tm_topo += "-".join([m4.group(1), m4.group(2)])
            dict_orf[orf.id] = orf


def write_results(output_file, format_long):
    global dict_orf
    
    out_handle = sys.stdout
    if output_file:
        out_handle = happyfile.hopen_write_or_else(output_file)
    
    if verbose:
        if output_file:
            xprint("Writing file: " + output_file)
        else:
            xprint("Writing stdout")

    if format_long:
        print >>out_handle, "\t".join(['id', 'family', 'subfamily', 'substrate', 'score', 'tc', 'tm_helices', 'tm_length', 'tm_topology'])
    else:
        print >>out_handle, "\t".join(['id', 'family', 'substrate', 'score'])
    
    for id in sorted(dict_orf.keys()):
        orf = dict_orf[id]
        
        orf.score = 0
        if orf.blast_cog:
            orf.score += weight_cog_positive
        if orf.blast_cog and orf.blast_cog in dict_fp_cog:
            orf.score += weight_cog_negative
        if orf.blast_tc:
            orf.score += weight_tc
        if orf.pfam_id:
            orf.score += weight_pfam
        if int(orf.tm_helices) >= min_tm_helices:
            orf.score += weight_tm_helices

        if orf.score >= min_orf_score:
            str_tc = str_family = str_subfamily = str_substrate = ""
            
            if orf.blast_cog in dict_cog_info:
                cog_info = dict_cog_info[orf.blast_cog]
                str_tc, str_family, str_subfamily, str_substrate = cog_info.tc, cog_info.family, cog_info.subfamily, cog_info.substrate
            elif orf.pfam_id in dict_hmm_info:
                hmm_info = dict_hmm_info[orf.pfam_id]
                str_tc, str_family, str_subfamily = hmm_info.tc, hmm_info.family, hmm_info.subfamily
            elif orf.blast_tc in dict_transporter:
                trans = dict_transporter[orf.blast_tc]
                str_tc, str_family, str_substrate = orf.blast_tc, trans.family, trans.substrate

            if str_family:
                if format_long:
                    print >>out_handle, "\t".join([orf.id, str_family, str_subfamily, str_substrate, str(orf.score), str_tc, str(orf.tm_helices), str(orf.tm_length), orf.tm_topo])
                else:
                    if str_subfamily and str_subfamily != str_family:
                        str_family = str_family + " (" + str_subfamily + ")"
                    print >>out_handle, "\t".join([orf.id, str_family, str_substrate, str(orf.score)])

def test_db():
    failed = 0
    table_found = []
    try:
        transdb_con = sqlite3.connect(trans_db_file)
        transdb_cursor = transdb_con.cursor()

        for row in transdb_cursor.execute('SELECT * from sqlite_master'):
            if row[0] == 'table':
                xprint("found table " + row[1])
                table_found.append(row[1])
    except sqlite3.Error:
        xprint("[test_db]: database error")
        failed = 1
    
    if len(set(table_found) & set(['cog_info', 'hmm_info', 'cog_seguid', 'fp_cog', 'tc_seguid', 'transporter'])) != 6:
        xprint("[test_db]: not all tables found in sqlite database")
        failed = 1

    if failed:
        xprint("[test_db]: failed")
    else:
        xprint("[test_db]: passed")

    return failed

def test_progs():
    failed = 0
    rc = subprocess.call('blastall &>/dev/null', shell=True)
    if rc == 1:
        xprint("blastall found.")
    else:
        xprint("[test_progs]: blastall not found.  Install or update path.")
        failed += 1
    
    rc = subprocess.call('hmmscan -h &>/dev/null', shell=True)
    if rc == 0:
        xprint("hmmscan found.")
    else:
        xprint("[test_progs]: hmmscan not found.  Install or update path.")
        failed += 1

    rc = subprocess.call('echo "X" | tmhmm &>/dev/null', shell=True)
    if rc == 0:
        xprint("tmhmm found.")
    else:
        xprint("[test_progs]: tmhmm not found.  Install or update path.")
        failed += 1

    if failed:
        xprint("[test_progs]: failed tests = " + str(failed))
    else:
        xprint("[test_progs]: passed")

    return failed

def test_all():
    failed = 0
    failed += test_db()
    failed += test_progs()

    if failed:
        xprint("[test_all]: failed tests = " + str(failed))
    else:
        xprint("[test_all]: all tests passed")

###

def main(argv):
    xprint("[find_transporters] start: " + time.asctime(time.localtime()))
    help = "\n".join([
        "find_transporters v0.1 (Apr 13, 2016)",
        "Membrane-bound transporter prediction",
        "",
        "Usage: "+argv[0]+" (options)",
        "   -f file        : FASTA file (required if missing any of -b,-p,-t)",
        "   -b file        : BLASTP output file (required if missing -f)",
        "   -p file        : HMM output file (required if missing -f)",
        "   -t file        : TMHMM output file (required if missing -f)",
        "",
        "   -c, --cpus int : parallel processes for BLASTP, HMMER, TMHMM (default: 3)",
        "   -h, --help     : help",
        "   -l             : extended output format",
        "   -m num         : minimum evidence score to output (default: 2)",
        "   -n             : only accept input files (-b,-p,-t); do NOT run processes",
        "   -o file        : output file (default: stdout)",
        "   -s             : brief output format (default)",
        "   -v, --verbose  : more information to stderr", ""])

    global cpus
    global min_orf_score
    global verbose
            
    fasta_file = ""
    blast_file = ""
    hmm_file = ""
    tmhmm_file = ""
    output_file = ""
    format_long = False
    do_not_run = False
    unused_args = []
    
    try:
        opts, args = getopt.getopt(argv[1:], "c:b:f:hlm:no:p:st:v", ["cpus=", "help", "verbose", "test"])
    except getopt.GetoptError:
        xprint(help)
        sys.exit(2)
    
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            xprint(help)
            sys.exit()
        elif opt == '--test':
            test_all()
            sys.exit()
        elif opt == '-f':
            fasta_file = arg
        elif opt == '-b':
            blast_file = arg
        elif opt == '-p':
            hmm_file = arg
        elif opt == '-t':
            tmhmm_file = arg
        elif opt in ("-c", "--cpus"):
            cpus = int(re.sub('=','', arg))
        elif opt == '-l':
            format_long = True
        elif opt == '-m':
            min_orf_score = float(re.sub('=','', arg))
        elif opt == '-n':
            do_not_run = True
        elif opt == '-o':
            output_file = arg
        elif opt == '-s':
            format_long = False
        elif opt in ("-v", "--verbose"):
            verbose = True
        else:
            unused_args.append(opt)

    if verbose:
        xprint("\n".join([
            "Program:        " + prog_path,
            "Database:       " + trans_db_file,
            "HMM models:     " + trans_hmm_file,
            "Blast DB:       " + trans_blast_file,
            "Score COG:      " + str(weight_cog_positive),
            "Score not COG:  " + str(weight_cog_negative),
            "Score tc:       " + str(weight_tc),
            "Score Pfam:     " + str(weight_pfam),
            "Score tm:       " + str(weight_tm_helices),
            "Min tm helices: " + str(min_tm_helices),
            "Max Blast e:    " + str(max_blast_evalue),
            "Max HMM e:      " + str(max_hmm_evalue),
            "CPUs:           " + str(cpus),
            "Min score:      " + str(min_orf_score), ""]))

    if len(unused_args) > 0:
        xprint(help)
        for arg in unused_args:
            xprint("Unused argument: " + str(arg))
        sys.exit()

    if do_not_run and not (blast_file and hmm_file and tmhmm_file):
        xprint(help + "\nBLAST, HMM, TMHMM files required with -n.")
        sys.exit()
            
    if not (fasta_file or (blast_file and hmm_file and tmhmm_file)):
        xprint(help + "\nFASTA file or output from all BLAST, HMM, TMHMM required.")
        sys.exit()

    if blast_file and hmm_file and tmhmm_file:
        if fasta_file:
            xprint("Ignoring FASTA file")
    else:
        # run blast, hmm, and/or tmhmm if outputs are not given
        if fasta_file:
            base_filename = os.path.basename(fasta_file)
            if not blast_file:
                blast_file = base_filename + ".blast"
                run_type_queue.append(Run_list_element(Run_type.blast, fasta_file, blast_file))
            if not hmm_file:
                hmm_file = base_filename + ".hmm"
                run_type_queue.append(Run_list_element(Run_type.hmm, fasta_file, hmm_file))
            if not tmhmm_file:
                tmhmm_file = base_filename + ".tmhmm"
                run_type_queue.append(Run_list_element(Run_type.tmhmm, fasta_file, tmhmm_file))
            if run_type_queue and not do_not_run:
                run_all_parallel()

    read_trans_db()
    read_blast(blast_file)
    read_hmms(hmm_file)
    read_tmhmm(tmhmm_file)
    
    write_results(output_file, format_long)
    xprint("[find_transporters] end: " + time.asctime(time.localtime()))

if __name__ == "__main__":
    main(sys.argv)
