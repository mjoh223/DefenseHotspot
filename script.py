from Bio import SeqIO
import glob
import pandas as pd
import os
from phamlite import drawOrf, graphing
from multiprocessing import Pool
import csv
import re
from Bio.SeqUtils import GC
import numpy as np
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
import ast
import subprocess

def findSystems(DHS, assembly, contig):
    mylist = list()
    print(len(DHS))
    for feature in DHS:
        id_='{}q{}_{}'.format(assembly, contig, feature.qualifiers['protein_id'][0].split('.')[0].replace('_',''))
        if id_ in dhs_family_dict:
            rep = dhs_family_dict[id_]
            mylist.append(rep)
        else:
            break
    print(mylist)
    return mylist

def createDHSprotfamily(filename):
    file1 = open(filename, 'r')
    lines = file1.readlines()
    mydict = dict()
    for line in lines:
        rep, mem = line.rstrip().split('\t')
        mydict[mem] = rep
    return mydict

def countSpaces(DHS):
    for i, gene in enumerate(DHS):
        try:
            print(DHS[i+1].location.start - gene.location.end)
            #print(np.absolute(gene.location.end - DHS[i+1].location.start))
        except:
            pass

def parseDomains(filename):
    file1 = open(filename, 'r')
    lines = file1.readlines()
    listOfLines = list()
    for i,line in enumerate(lines):
        if line[0] != '#':
            items = line.strip()
            items = items.split('\t')
            listOfLines.append(items)
    di = {}
    for i, line in enumerate(listOfLines):
        if line[0] == 'QUERY':
            query = line[4]
            #query = '|'.join(query.split(',')[1:])
        if line == ['DOMAINS']:
            c = 1
            while listOfLines[i+c] != ['ENDDOMAINS']:
                hit = listOfLines[i+c]
                hit = '|'.join([hit[6],hit[8],hit[9]])
                if query in di.keys():
                    di[query][0].add(hit)
                else:
                    di[query] = [set([hit])]#[hit[9],]
                c +=1
    for k,v in di.items():
        v = list(v[0])
        v = [x.split('|') for x in v]
        [print('{}\t{}\t{}\t{}'.format(k,*x)) for x in v]
        #domains = ['{}\t{}\t{}'.format(x) for x in list(v)[0].split('|')]
        #print('{}\t'+domains)
    return di

def domainAnalysis(filename):
    command = ['rpsblast '+ '-query '+ filename+' '+ '-db ' '/hdd-roo/DHS/pfam/db/Pfam '+ '-evalue '+ '0.01 '+ '-outfmt '+ '11 '+ '-out '+ 'DHS_domains.txt ']
    subprocess.run(command, stdout=subprocess.PIPE)

def gc_content(seq):
    return round((seq.count('C') + seq.count('G')) / len(seq) * 100)

def gc_content_subsec(seq, k=1000):
    res = []
    for i in range(0, len(seq) - k+1, k):
        subseq = seq[i:i+k]
        res.append(gc_content(subseq))
    return res

def plotGC(filename):
    file_ = open(filename, 'r')
    lines = file_.readlines()
    list_of_gc = [x.rstrip().strip('][').split(', ') for x in lines]
    for gc in list_of_gc:
        gc = gc[:-10]
        if len(gc) < 400:
            mylist = list()
            for i in gc:
                try:
                    mylist.append(int(i))
                except:
                    1+1
            y = np.array(mylist)
            #x = np.array(range(len(mylist)))
            plt.plot(y,alpha=0.2,color='green')
    plt.ylim([30, 80])
    plt.savefig("gc.png")


def regressionGC(seq):
    y = np.array(gc_content_subsec(seq[1000:-10000]))
    #print('{}'.format([x for x in y]))
    x = np.array(range(len(y))).reshape((-1,1))
    model = LinearRegression()
    model.fit(x,y)
    r_sq = model.score(x, y)
    print('coefficient of determination:', r_sq)
    print('slope:', model.coef_)

def calculateGC(record, DHS, a ,b ,c):
    fna = str(record.seq)
    first = int(DHS[0].location.start)
    last = int(DHS[-1].location.end)
    print('{}\t{}\t{}\t{}\t{}'.format(a,b,c,GC(fna),GC(fna[first:last])))

def parallelcreategembasefromDHS(DHS_record_features, assembly, contig):
    f = open(os.path.join('/hdd-roo/DHS/DHS_proteins/', contig), "w")
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

def graph(systems, DHS, contig):
    traces = drawOrf(systems, DHS)
    first = int(DHS[0].location.start)
    last = int(DHS[-1].location.end)
    size = last - first
    graphing(traces, contig, size)

def grabDHS(record, boundary_set):
    l_coord = [x[0] for x in enumerate(record.features) if x[1].qualifiers.get('locus_tag', [''])[0] in boundary_set[0]]
    r_coord = [x[0] for x in enumerate(record.features) if x[1].qualifiers.get('locus_tag', [''])[0] in boundary_set[1]]
    if len(l_coord) == 2 and len(r_coord) == 2:
        assembly = record.dbxrefs[0].split(':')[1].replace('_','').replace('.','')
        contig = record.name.replace('_','').replace('.','')
        coords = sorted([l_coord,r_coord])
        offset = 10 # number of features to extend the boundary for each side
        DHS = record.features[coords[0][0]-offset:coords[1][1]+offset]
        desc = record.description
        DHS = [x for x in DHS if x.type in ['CDS', 'tRNA']] #and 'protein_id' in x.qualifiers]
        if len(DHS) < 148+11:#'aeruginosa' in desc and len(DHS) > 0:
            #first = int(DHS[0].location.start)
            #last = int(DHS[-1].location.end)
            #fna = str(record.seq)[first:last]
            if 'tRNA' not in [x.type for x in DHS][:int(len(DHS)/2)]: #make sure the tRNA side is the beginning on the list
                DHS.reverse()
            #    fna = str(record.seq[first:last].reverse_complement())
            #countSpaces(DHS)
            #DHS = [x for x in DHS if x.qualifiers['protein_id'][0].replace('_','').split('.')[0] not in systems]
            #reps = findSystems(DHS, assembly, contig)
            #parallelcreategembasefromDHS(DHS, assembly, contig)
            #regressionGC(fna)
            #first = int(DHS[0].location.start)
            #last = int(DHS[-1].location.end)
            #print('{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(desc, assembly, contig, first, last, last-first,len(DHS)))
            graph(systems, DHS, contig)
        #f = open(os.path.join('/hdd-roo/DHS/parallel/', contig), "w")
        #dhs_len = len(list(filter(lambda x: x.type == 'CDS', DHS)))
        #first = DHS[0].qualifiers['locus_tag'][0]
        #last = DHS[-1].qualifiers['locus_tag'][0]
        #f.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(desc, assembly, contig, dhs_len, first, last))
        #f.close()
            #parallelcreategembasefromDHS(DHS[1:-1], assembly, contig)
        #graph(systems, DHS, contig)
        #print(lol)
def readGenbanks(filename, boundary_set):
    gb_files = glob.glob(filename)
    for gb_file in gb_files:
        [grabDHS(gb_record, boundary_set) for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank")]

def readISLANDproteins(filename):
    reader = list(csv.reader(open(filename)))
    for row in reader:
        systems = row
        break
    mydict = dict()
    for row in reader[1:]:
        for i, cell in enumerate(row[2:]):
            prots = re.split(':|&', cell) 
            for prot in prots:
                if prot is not '':
                    mydict[prot] = systems[i+2]
    return mydict

def readMMseqOut(filenames):
    df0 = pd.read_csv(filenames[0], delimiter = '\t', header=None)
    df1 = pd.read_csv(filenames[1], delimiter = '\t', header=None)
    df0 = df0[(df0.iloc[:,2] > 0.8) & (df0.iloc[:,3] > 93.6)] 
    df1 = df1[(df1.iloc[:,2] > 0.8) & (df1.iloc[:,3] > 1678.4)] 
    df0 = set(df0.iloc[:,1])
    df1 = set(df1.iloc[:,1])
    return [df0, df1]

def runner(gb_file):
    [grabDHS(gb_record, boundary_set) for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank")]

def parallel(genbanks):
    #for genbank in genbanks:
    p = Pool(128)
    p.map(runner, genbanks)

boundary_set = readMMseqOut(['/hdd-roo/DHS/lb_results', '/hdd-roo/DHS/rb_results'])
genbanks = glob.glob('/hdd-roo/DHS/gbk/*gbk')
systems = readISLANDproteins('/hdd-roo/DHS/042822_island_results_prot.csv')
#parseDomains('/hdd-roo/DHS/pfam/050322_domain_results.tsv')
dhs_family_dict = createDHSprotfamily('/hdd-roo/DHS/predict_systems/db_cluster.tsv')
parallel(genbanks)
#plotGC('GC_values.txt')
#createGembase(glob.glob('/hdd-roo/DHS/gbk/*gbk'))
#readGenbanks('/hdd-roo/DHS/gbk/*gbk', boundary_set)
