#! /usr/bin/env python
import os
import re
import sys
import string
import optparse
scr_dir = os.path.abspath(os.path.dirname(sys.argv[0]))
sys.path.append(scr_dir)
import fasta
import math

def format(seq, N=60):
    nseg = int(math.ceil(len(seq)/(N+0.0)))
    return '\n'.join([seq[i*N:(i+1)*N] for i in range(nseg)])

parser = optparse.OptionParser(version="%prog ")

parser.add_option("-i", "--input", dest="input_fasta",action="store",
        help="name of input file in fasta format")
parser.add_option("-L", "--LibHmm", dest="hmm_path",action="store",
        default="HMMs",help="path of hmm database")
parser.add_option("-o", "--output", dest="out_fname",action="store",
        help="name of output file")
parser.add_option("-k", "--kingdoms", dest="kingdoms",action="store",
        default="arc,bac",help="kingdom used")
parser.add_option("-m", "--moltypes", dest="moltypes",action="store",
        default="lsu,ssu,tsu",help="molecule type detected")
parser.add_option("-e","--Evalue", dest="evalue",action="store",type="float",
        default=0.01,help="evalue cut-off for hmmsearch")
parser.add_option("-H","--hmmerhead",dest="hmmerhead",action="store",
        default="no",help="whether use HMMERHEAD to speedup hmmsearch")

try:
    (options, args) = parser.parse_args()
except:
    parser.print_help()
    sys.exit(1)

if options.input_fasta is None or options.hmm_path is None:
    parser.print_help()
    sys.exit(1)

try:
    os.environ["HMMERDB"] += ":"+os.path.abspath(options.hmm_path)
except: # because the machine may not have HMMERDB variable
    os.environ["HMMERDB"] = os.path.abspath(options.hmm_path)
#print os.environ["HMMERDB"]
out_fname = os.path.abspath(options.out_fname)
out_dir=os.path.dirname(out_fname)
fname = os.path.abspath(options.input_fasta)

tr = string.maketrans("gatcryswkmbdhvnGATCRYSWKMBDHVN","ctagyrswmkvhdbnCTAGYRSWMKVHDBN")


def rev_record(record):
    return ">"+record.header+"|rev\n"+format(record.sequence[::-1].translate(tr))

def parse_hmmsearch(src):
# function to parse hmmsearch output
    resu = []
    re_header = re.compile(r"^Query(\s+(sequence|HMM))?:\s+")
    re_domain = re.compile(r"^Parsed for domains:")
    re_model = re.compile(r"^Model|Sequence\s+Domain")
    re_match = re.compile(r"(\S+)\s+(\d+)\/(\d+)\s+(\d+)\s+(\d+).+?(\d+)\s+(\d+)\s+\S+\s+(\S+)\s+(\S+)\s*$")
    re_begin = re.compile(r"^\-\-")
    re_end = re.compile(r"^Alignments of")
    ready = -1
    if type(src) == str: # provide a filename:
        src = open(src)#.readlines()
    for line in src:
        if re_header.search(line):
            Query = line.strip().split()[-1]
            continue
        elif re_domain.search(line):
            ready = 0
            continue
        if re_model.search(line):
            ready = 1
        elif ready > 0 and re_begin.search(line):
            ready = 2
        elif re_end.search(line):
            break
        elif ready == 2:
            if re_match.search(line):
                vals = re_match.search(line).groups()
                resu.append('\t'.join([Query,vals[3],vals[4],vals[0],vals[5],vals[6],vals[7],vals[8],vals[0]]))
    return resu
    
records = [rec for rec in fasta.fasta_itr(fname)]
headers = [[rec.header,len(rec.sequence)] for rec in records]

ff = open(out_fname+'.fa','w')
for (i, rec) in enumerate(records):
    ff.write('>s'+str(i+1)+'\n'+format(rec.sequence)+'\n')
    ff.write('>s'+str(i+1)+'|rev\n'+format(rec.sequence[::-1].translate(tr))+'\n')
ff.close()
# a temporary fasta file, use s(int) to easy the parsing

hmm_resu = []
for kingdom in options.kingdoms.split(','):
    for moltype in options.moltypes.split(','):
#        print kingdom, moltype
        if options.hmmerhead == "no":
            cmd = 'hmmsearch -E %g %s_%s.fs.hmm %s' % (options.evalue,kingdom,moltype,out_fname+'.fa')
        else:
            cmd = 'hmmsearch --hmmerhead -E %g %s_%s.fs.hmm %s' % (options.evalue,kingdom,moltype,out_fname+'.fa')
#        print cmd
        hmm_resu += parse_hmmsearch(os.popen(cmd))

os.remove(out_fname+'.fa')

ff = open(out_fname,"w")
dict_rRNA = {'arc_lsu':'23S_rRNA','arc_ssu':'16S_rRNA','arc_tsu':'5S_rRNA',
             'bac_lsu':'23S_rRNA','bac_ssu':'16S_rRNA','bac_tsu':'5S_rRNA',
             'euk_lsu':'28S_rRNA','euk_ssu':'18S_rRNA','euk_tsu':'8S_rRNA'}

dict_read2kingdom = {}
for line in hmm_resu:
    [feature_type, r_start, r_end, read, h_start, h_end, score, evalue, read1] = line.strip().split('\t')
    read = read.split('|')[0]
    evalue = string.atof(evalue)
    kingdom = feature_type.split('_')[0]
    if read in dict_read2kingdom:
        if evalue < dict_read2kingdom[read][1]:
            dict_read2kingdom[read] = [kingdom, evalue]
    else:
        dict_read2kingdom[read] = [kingdom, evalue]

header = ['##seq_name','method','feature','start','end','evalue','strand','frame','attribute']
ff.write('\t'.join(header)+'\n')
for line in hmm_resu:
    [feature_type, r_start, r_end, read, h_start, h_end, score, evalue, read1] = line.strip().split('\t')
#    if string.atof(evalue) < evalue_cut and feature_type.split('_')[0] == orgn:
    if dict_read2kingdom[read.split('|')[0]][0] != feature_type.split('_')[0]:
        continue
    feature_type = dict_rRNA[feature_type]
    if read.endswith('|rev'):
        strand = '-'
        tmp = map(string.atoi,[r_start,r_end])
        pos = string.atoi(read[1:-4])-1
        header = headers[pos][0]
        L = headers[pos][1]
        [r_end,r_start] = [str(L+1-x) for x in tmp]
    else:
        strand = '+'
        pos = string.atoi(read[1:])-1
        header = headers[pos][0]
    ff.write('\t'.join([header, 'rna_hmm_fs','rRNA',r_start,r_end,evalue,strand,'NA',feature_type])+'\n')
ff.close()

