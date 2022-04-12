from Bio import SeqIO
import glob
import pandas as pd
import os
from phamlite import drawOrf, graphing
from multiprocessing import Pool
def parallelcreategembasefromDHS(DHS_record_features, assembly, contig):
    f = open(os.path.join('/hdd-roo/DHS/parallel/', contig), "w")
    for feature in DHS_record_features:
        if feature.type == 'CDS' and 'protein_id' in feature.qualifiers:
            f.write('>{}q{}_{}\n'.format(
                        assembly,
                        contig,
                        feature.qualifiers['protein_id'][0].split('.')[0].replace('_','')
                        ))
            f.write('{}\n'.format(feature.qualifiers['translation'][0]))

def createGembasefromDHS(DHS_record_features, assembly, contig):
    for feature in DHS_record_features:
        if feature.type == 'CDS' and 'protein_id' in feature.qualifiers:
            print('>{}q{}_{}'.format(
                        assembly,
                        contig,
                        feature.qualifiers['protein_id'][0].split('.')[0].replace('_','')
                        ))
            print(feature.qualifiers['translation'][0])

def createGembase(gb_files):
    for gb_file in gb_files:
        for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
            assembly = gb_record.dbxrefs[0].split(':')[1].replace('_','').replace('.','')
            for feature in gb_record.features:
                if feature.type == 'CDS' and 'protein_id' in feature.qualifiers:
                    print('>{}q{}_{}'.format(
                        assembly,
                        gb_record.name.replace('_','').replace('.',''),
                        feature.qualifiers['protein_id'][0].split('.')[0].replace('_','')
                        ))
                    print(feature.qualifiers['translation'][0])

def graph(DHS):
    traces = drawOrf(DHS)
    graphing(traces)

def grabDHS(record, boundary_set):
    l_coord = [x[0] for x in enumerate(record.features) if x[1].qualifiers.get('locus_tag', [''])[0] in boundary_set[0]]
    r_coord = [x[0] for x in enumerate(record.features) if x[1].qualifiers.get('locus_tag', [''])[0] in boundary_set[1]]
    if len(l_coord) == 2 and len(r_coord) == 2:
        assembly = record.dbxrefs[0].split(':')[1].replace('_','').replace('.','')
        contig = record.name.replace('_','').replace('.','')
        coords = sorted([l_coord,r_coord])
        DHS = record.features[coords[0][0]:coords[1][1]]
        #print(DHS)
        parallelcreategembasefromDHS(DHS, assembly, contig)
        #graph(DHS)
        #print(lol)
def readGenbanks(filename, boundary_set):
    gb_files = glob.glob(filename)
    for gb_file in gb_files:
        [grabDHS(gb_record, boundary_set) for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank")]

def readMMseqOut(filenames):
    df0 = pd.read_csv(filenames[0], delimiter = '\t', header=None)
    df1 = pd.read_csv(filenames[1], delimiter = '\t', header=None)
    df0 = set(df0.iloc[:,1])
    df1 = set(df1.iloc[:,1])
    return [df0, df1]

def runner(gb_file):
    [grabDHS(gb_record, boundary_set) for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank")]

def parallel(genbanks):
    #for genbank in genbanks:
    p = Pool(127)
    p.map(runner, genbanks)

boundary_set = readMMseqOut(['/hdd-roo/DHS/lb_results', '/hdd-roo/DHS/rb_results'])
genbanks = glob.glob('/hdd-roo/DHS/gbk/*gbk')
parallel(genbanks)
#createGembase(glob.glob('/hdd-roo/DHS/gbk/*gbk'))
#readGenbanks('/hdd-roo/DHS/gbk/*gbk', boundary_set)


