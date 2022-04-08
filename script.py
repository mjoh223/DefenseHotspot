from Bio import SeqIO
import glob
import pandas as pd

def readGenbanks(filename):
    return SeqIO.parse(open(filename,"r"), "genbank")

def readMMseqOut(filenames):
    df0 = pd.read_csv(filenames[0], delimiter = '\t', header=None)
    df1 = pd.read_csv(filenames[1], delimiter = '\t', header=None)
    df0 = set(df0.iloc[:,1])
    df1 = set(df1.iloc[:,1])
    return df0 | df1
boundary_set = readMMseqOut(['/hdd-roo/DHS/lb_results', '/hdd-roo/DHS/rb_results'])
print(boundary_set)
