from Bio import SeqIO
import glob
import pandas as pd


def grabDHS(record, boundary_set):
    l_coord = [x[0] for x in enumerate(record.features) if x[1].qualifiers.get('locus_tag', [''])[0] in boundary_set[0]]
    r_coord = [x[0] for x in enumerate(record.features) if x[1].qualifiers.get('locus_tag', [''])[0] in boundary_set[1]]
    if len(l_coord) > 0 and len(r_coord) > 0:
        print([
        print(record.features[l_coord[0]:r_coord[0]])

def readGenbanks(filename, boundary_set):
    gb_files = glob.glob(filename)
    for gb_file in gb_files:
        DHS = [grabDHS(gb_record, boundary_set) for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank")]

def readMMseqOut(filenames):
    df0 = pd.read_csv(filenames[0], delimiter = '\t', header=None)
    df1 = pd.read_csv(filenames[1], delimiter = '\t', header=None)
    df0 = set(df0.iloc[:,1])
    df1 = set(df1.iloc[:,1])
    return [df0, df1]
boundary_set = readMMseqOut(['/hdd-roo/DHS/lb_results', '/hdd-roo/DHS/rb_results'])
readGenbanks('/hdd-roo/DHS/gbk/*gbk', boundary_set)
