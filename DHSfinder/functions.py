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
# 1- preprocessing: split gembase into 128 parts
# 2- functions are ran on each gembase partition

def preprocess(user_gembase, working_dir):
    with open(user_gembase, 'a') as outFile:
    record_ids = set()
    for record in SeqIO.parse('input.fasta', 'fasta'):
        if record.id not in record_ids:
            record_ids.add(record.id)
            SeqIO.write(record, outFile, 'fasta')
            
if __name__ == "__main__":
    working_dir = '/hdd-roo/DHS/DHSfinder_out/'
    preprocess('/hdd-roo/DHS/short_gembase.fa', working_dir)
print(i)

def calc_batch_sizes(n_tasks: int, n_workers: int) -> list:
    x = int(n_tasks / n_workers)
    y = n_tasks % n_workers
    batch_sizes = [x + (y > 0)] * y + [x] * (n_workers - y)
    return batch_sizes

def build_batch_ranges(batch_sizes: list) -> list:
    upper_bounds = [*itertools.accumulate(batch_sizes)]
    lower_bounds = [0] + upper_bounds[:-1]
    batch_ranges = [range(l, u) for l, u in zip(lower_bounds, upper_bounds)]
    return batch_ranges

def cluster(gem_file):
    #input: gembase file
    #returns: classes with gene families
    working_path = '/hdd-roo/DHS/DHSfinder_out/'
    binb = '/home/jbd_apps/dash_app/MMseqs2/build/bin'
    cluster = ['{}/mmseqs'.format(binb),
                'easy-cluster',
                gem_file,
                os.path.join(working_path, 'clu'),
                os.path.join(working_path, 'tmp')]
    #subprocess.run(cluster, shell=False)
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
    gem_file = '/hdd-roo/DHS/pae_short_gem.fa'#'/hdd-roo/DHS/p_aeruginosa_total_gembase.txt'
    print('reading cluster tsv')
    tsv = cluster(gem_file)
    print('making lookup table')
    cluster_dict = get_cluster_dict(tsv) #returns a dictionary of the mmseqs cluster representatives
    print('getting contig-protein dict')
    parallel_contig_dict(gem_file, cluster_dict)
    #contig_protein_dict, db_size = get_contig_dict(gem_file, cluster_dict)
    print('done')
    #contig_list = contig_protein_dict.items()
    #regions_dict = parallel_regions(contig_list)
    #tophits(regions_dict)
