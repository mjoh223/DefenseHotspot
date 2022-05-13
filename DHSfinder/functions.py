from classes import *
from multiprocessing import Pool
import itertools
import csv
import pandas as pd
import multiprocessing
from Bio import SeqIO
import numpy as np
import os
import pickle
import subprocess
import uuid
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import pyfasta

def parseGemVfile(gem_file):
    #return a dictionary contig for each assembly in gembase
    assemblies = set()
    contig_protein_dict = {}
    index = 0
    all_contigs = list()

def cluster(gem_file):
    #input: gembase file
    #returns: classes with gene families
    working_path = '/hdd-roo/DHS/scripts/DHSfinder_out/'
    binb = '/home/jbd_apps/dash_app/MMseqs2/build/bin'
    cluster = ['{}/mmseqs'.format(binb),
                'easy-cluster',
                gem_file,
                os.path.join(working_path, 'clu'),
                os.path.join(working_path, 'tmp')]
    subprocess.run(cluster, shell=False)
    return pd.read_csv(os.path.join(working_path, 'clu_cluster.tsv'), sep='\t',header=None, names=['representative','member'])

def get_cluster_dict(clust_pd):
    return dict(zip(clust_pd.member, clust_pd.representative))

def get_contig_dict(filename, cluster_dict):
    contig_protein_dict = {}
    handle = open(filename, 'r')
    index = 0
    current_contig = ''
    genome_set = set()
    for record in SeqIO.parse(handle, 'fasta'):
        genome_set.add(record.id.split('q')[0])
        contig, identifier = record.id.split('q')[1].split('_')
        if contig not in contig_protein_dict:
            contig_protein_dict[contig] = []
        contig_protein_dict[contig].append(SeqUnit(contig=contig, position=index, sequence_type='prot', sequence=record.seq, family=cluster_dict[record.id], identifier=identifier))
        index += 1
    return contig_protein_dict, len(genome_set)

def getregionsparallel(contig_tuple):
    regions_dict = {}
    win = range(3,10)
    contig_id, contig = contig_tuple
    for i, gene in enumerate(contig):
        for reg in win:
            if i > len(contig)-reg:
                break
            if len(contig) <= i+reg - i:
                break
            potentials = contig[i:i+reg]
            rid = RegionID(potentials[0].family, potentials[-1].family)
            if rid in regions_dict:
                regions_dict[rid].regions.append(potentials[1:-1])
            else:
                regions_dict[rid] = AcessoryRegion(potentials[0], potentials[-1])
                regions_dict[rid].regions.append(potentials[1:-1])
    return regions_dict

def tophits(regions_dict):
    hits = dict()
    for k,v in regions_dict.items():
        myset = set()
        for accreg in v:
            for regs in accreg.regions:
                for gene in regs:
                    myset.add(gene.family)
        hits[k] = myset
    sorted_lis = list(sorted(hits.items(), key=lambda item: len(item[1]), reverse=True))
    for info in sorted_lis[:5]:
        rid = info[0]
        acc = info[1]
        print(rid.first_bound)
        print(rid.second_bound)
        print(acc)

def parallel_contig_dict(gem_file, cluster_dict):
    p = Pool(128)
    p.map()

def parallel_regions(contig_list):
    print(len(contig_list))
    p = Pool(128)
    print('parallelizing boundary search')
    dicts = p.map(getregionsparallel, contig_list)
    print('merging')
    dd = defaultdict(list)
    for d in dicts:
        for key, value in d.items():
            dd[key].append(value)
    print(len(dd))
    return dd

if __name__ == "__main__":
    gem_file = '/hdd-roo/DHS/short_gembase.fa'#'/hdd-roo/DHS/p_aeruginosa_total_gembase.txt'
    tsv = cluster(gem_file)
    print('filtering potential DHS boundaries')
    cluster_dict = get_cluster_dict(tsv) #returns a dictionary of the mmseqs cluster representatives
    print('getting contig protein dict')
    contig_protein_dict, db_size = get_contig_dict(gem_file, cluster_dict)
    contig_list = contig_protein_dict.items()
    regions_dict = parallel_regions(contig_list)
    tophits(regions_dict)
