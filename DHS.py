import itertools
import csv
import pandas as pd
import multiprocessing
from Bio import SeqIO
from scipy import stats
import numpy as np
from matplotlib import colors
import os
import pickle
import subprocess
import seaborn as sns
from matplotlib import pyplot as plt
import uuid
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import islice
from Bio import Entrez
import random
import drawSvg as draw
# colors_dict = colors.CSS4_COLORS
#
#this defines the SeqUnit and AcessoryRegion classes


class SeqUnit:

    def __init__(self, contig, sequence, position, sequence_type, family, identifier):

        self.sequence = sequence
        self.position = int(position)
        self.contig = contig
        self.sequence_type = sequence_type
        self.family = family
        self.identifier = identifier

    # if the sequence type is nuc then the position is the actual position in the genome of the kmer and if the sequence type is protein then the distance
    # is the number of genes

    def can_it_join(self, other_sequnit, min_space, max_distance):

        return (self.contig == other_sequnit.contig and abs(self.position - other_sequnit.position) < max_distance and abs(self.position - other_sequnit.position) > min_space)


class AcessoryRegion:

    def __init__(self, first_bound_family, second_bound_family):

        self.first_bound_family = first_bound_family
        self.second_bound_family = second_bound_family
        self.regions = []
        self.score = 0


class RegionID:

    def __init__(self, first_bound, second_bound):

        self.first_bound = first_bound
        self.second_bound = second_bound

    def __eq__(self, other):

        return (self.first_bound == other.first_bound or self.first_bound == other.second_bound) and (self.second_bound == other.first_bound or self.second_bound == other.second_bound)

    def __hash__(self):

        return hash(self.first_bound) ^ hash(self.second_bound)

def region_conservation(regions):

    #regions should be a list of lists and this will return the average similarity of every pair of lists in the list

    x = list(itertools.chain.from_iterable(regions))
    y = [thing.family for thing in x]
    return len(set(y)) / len(y)

def extract_sequences(genome, bound1, bound2, dictionary):

    region = []

    index1 = bound1.split('_')[1]
    index2 = bound2.split('_')[2]

    for thing in dictionary[genome][index1:index2 + 1]:

        region.append(thing)

    return region

def get_cluster_dict():
    filename = '/hdd-roo/genomes_master/contig_gembases/Pseudomonas_aeruginosa/p_aeruginosa_total_gembase_clust25_cluster.tsv'#'/hdd-roo/genomes_master/contig_gembases/paer/paer_total_gembase_clust_cluster.tsv'
    handle = open(filename)
    read_tsv = csv.reader(handle, delimiter="\t")
    cluster_dict = {}
    for rep, mem in read_tsv:
        cluster_dict[mem] = rep
    handle.close()
    return cluster_dict

def get_contig_dict(cluster_dict):
    contig_protein_dict = {}
    handle = open('/hdd-roo/genomes_master/contig_gembases/Pseudomonas_aeruginosa/p_aeruginosa_total_gembase.txt', 'r')
    index = 0
    current_contig = ''
    genome_set = set()
    for record in SeqIO.parse(handle, 'fasta'):
        if 'q' in record.id:
            genome_set.add(record.id.split('q')[0])
            contig = record.id.split('q')[1].split('_')[0]
            identifier = record.id.split('q')[1].split('_')[1]
            if contig not in contig_protein_dict:
                contig_protein_dict[contig] = []
            contig_protein_dict[contig].append(SeqUnit(contig=contig, position=index, sequence_type='prot', sequence=record.seq, family=cluster_dict[record.id], identifier=identifier))
            index += 1
    print(genome_set)
    return contig_protein_dict, len(genome_set)

def get_regions2(contig_protein_dict):
    n = 0
    regions_dict = {}
    win = range(3,100)
    for contig_id, contig in contig_protein_dict.items():
        #if n % 500 == 0:
        #    print(n)
        if n == 10000:
            break
        for i, gene in enumerate(contig):
            for reg in win:
                if i > len(contig)-reg:
                    break
                if len(contig) <= i+reg - i:
                    break
                potentials = contig[i:i+reg]
                rid = RegionID(potentials[0].family, potentials[-1].family)
                #rid = tuple(sorted([potentials[0].family, potentials[-1].family]))
                if rid in regions_dict:
                    regions_dict[rid].regions.append(potentials)
                else:
                    regions_dict[rid] = AcessoryRegion(potentials[0], potentials[-1])
                    regions_dict[rid].regions.append(potentials)
        n+=1
    return regions_dict

def get_regions(contig_protein_dict):
    regions_dict = {}
    max_distance = 10
    print('getting regions...')
    index_big = 0
    print(len(contig_protein_dict))
    for contig_id, contig in contig_protein_dict.items():
        if index_big % 100000 == 0:
            print(index_big)
        n = len(contig)
        print(n)
        for i, unit in enumerate(contig): #unit is a gene
            lower_bound = i #- max_distance
            upper_bound = i + max_distance
            potentials = []
            if n < 2 * max_distance:
                potentials = contig
            elif upper_bound > n - 1:
                potentials = contig[lower_bound:]
            #this allows the region to be 2x max_distance, which is fine just make sure you are careful
            else:
                potentials = contig[lower_bound:upper_bound]
            for j in range(i,len(potentials)):
                current_id = RegionID(unit.family, potentials[j].family)
                if current_id in regions_dict:
                    regions_dict[current_id].regions.append(potentials)
                else:
                    ar = AcessoryRegion(unit.family, potentials[j].family)
                    ar.regions = potentials
                    regions_dict[current_id] = ar
        index_big += 1
        break
    return regions_dict

def eric():
    regions_dict = {}
    max_distance = 10
    print('here')
    index_big = 0
    print(len(contig_protein_dict))
    for contig in contig_protein_dict:
        if index_big % 500 == 0:
            print(index_big)
        for i in range(len(contig_protein_dict[contig])):
            lower_bound = i - max_distance
            upper_bound = i + max_distance
            potentials = []
            if len(contig_protein_dict[contig]) < 2 * max_distance:
                potentials = contig_protein_dict[contig]
            elif lower_bound < 0:
                potentials = contig_protein_dict[contig][0:upper_bound + 1]
            elif upper_bound > len(contig_protein_dict[contig]) - 1:
                potentials = contig_protein_dict[contig][lower_bound:]
            #this allows the region to be 2x max_distance, which is fine just make sure you are careful
            else:
                potentials = contig_protein_dict[contig][lower_bound:upper_bound]
            current_protein = contig_protein_dict[contig][i]
            index = potentials.index(current_protein)
            for j in range(len(potentials)):
                protein2 = potentials[j]
                if protein2 != current_protein:
                    current_id = RegionID(current_protein.family, protein2.family)
                    if current_id not in regions_dict:
                        regions_dict[current_id] = AcessoryRegion(current_protein.family, protein2.family)
                    if j > index:
                        regions_dict[current_id].regions.append(potentials[index:j] + [potentials[j]])
                    else:
                        regions_dict[current_id].regions.append(potentials[j:index] + [potentials[index]])
        index_big += 1

def read_ISLAND():
    from operator import methodcaller
    import numpy as np
    filename = '/hdd-rabbit/jbd_apps/final_ISLAND_proteins_4_7_FZ_and_smaller_window.csv'
    with open(filename) as csvfile:
        rows = csv.reader(csvfile, delimiter=',')
        hits = np.array(list(rows))
    arr = hits[:,3:].transpose()
    mydict = {}
    for x in arr:
        for y in x[1:]:
            if y is not '':
                pp = y.split(':')
                for p in pp:
                    mydict[p] = x[0].rstrip()
    return mydict
    #print(list(map(methodcaller("split", ":"), arr)))

def DHS(regions_dict, island_dict, db_size, C_thresh):
    #print('systems\tconservation\thyper\twidth\t,width_var\tperc_system')
    hs = {}
    for k,v in regions_dict.items():
        #if k.first_bound == 'GCA001879525qCP013113.1_APC72375.1' or k.second_bound == 'GCA001879525qCP013113.1_APC72375.1' and k.first_bound == 'GCA900618265qLR130531.1_VDL17998.1' or k.second_bound == 'GCA900618265qLR130531.1_VDL17998.1':
        #    desc = 'DHS1'
        #else:
        #    desc = 'notDHS1'
        #Cscore = len(set([x[0].contig for x in v.regions]))/db_size
        #print(Cscore)
        #mean_width = sum([len(x) for x in v.regions]) / len(v.regions)
        #width_variance = sum([(len(x)-mean_width)**2 for x in v.regions])/len(v.regions)
        if True: #width_variance < 1 and Cscore >= 0.8:
            systems = set()
            diversity = set()
            n,s = 0, 0  
            for region in v.regions:
                n += len(region)
                diversity.update([x.family for x in region])
                if any([gene.identifier in island_dict for gene in region]): #how many regions have atleast one system
                    s += 1 
                for gene in region:
                    if gene.identifier in island_dict:
                        systems.add(island_dict[gene.identifier])
            Hscore = len(diversity)/n
            perc_system = s/len(v.regions)
            if len(systems) > 5:
                hs[k] = v
                print(systems)
                #print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(len(systems), Cscore, Hscore, mean_width, width_variance, perc_system, desc))
    return hs
            #print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
            #    k.first_bound, 
            #    k.second_bound, 
            #    len(systems), 
            #    Cscore, 
            #    Hscore, 
            #    mean_width, 
            #    width_variance, 
            #    perc_system,
            #    str(systems),
            #    str(diversity)))

def DHS2(hs):
    for k,v in hs.items():
        a=k.first_bound.split('_')[1]
        b=k.second_bound.split('_')[1]
        if k.first_bound == 'GCA001879525qCP013113.1_APC72375.1' or k.second_bound == 'GCA001879525qCP013113.1_APC72375.1' and k.first_bound == 'GCA900618265qLR130531.1_VDL17998.1' or k.second_bound == 'GCA900618265qLR130531.1_VDL17998.1':
            desc = 'DHS1'
        else:
            desc = 'notDHS1'
        Cscore = len(set([x[0].contig for x in v.regions]))/db_size
        mean_width = sum([len(x) for x in v.regions]) / len(v.regions)
        width_variance = sum([(len(x)-mean_width)**2 for x in v.regions])/len(v.regions)
        systems = set()
        diversity = set()
        n,s = 0, 0
        filename = '/hdd-rabbit/jbd_apps/candidte_DHS_unq/' + "_".join([desc,a,b]) + '.txt'
        #f = open(filename, "w")
        for region in v.regions:
            n += len(region)
            diversity.update([x.family for x in region])
            if any([gene.identifier in island_dict for gene in region]): #how many regions have atleast one system
                s += 1
            print([domain_dict.get(x.family,'no') for x in region])
            for gene in region: 
                if gene.identifier in island_dict:
                    systems.add(island_dict[gene.identifier])
        #for family in diversity:
        #    for region in v.regions:
        #        hit = False
        #        for gene in region:
        #            if family == gene.family:
        #                f.write('>{}\n'.format(gene.family))
        #                f.write(str(gene.sequence))
        #                f.write('\n')
        #                hit = True
        #                break
        #        if hit:
        #            break
        #f.close()
        

        
        try:
            ainfo = getinfo(a)[0].get('Title', 'nf')
        except:
            ainfo = 'nf'
        try:
            binfo = getinfo(b)[0].get('Title', 'nf')
        except:
            binfo = 'nf'

        Hscore = len(diversity)/n
        perc_system = s/len(v.regions)
        #print([k.first_bound, k.second_bound])
        #print(systems)
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t".format(len(systems), Cscore, Hscore, mean_width, width_variance, perc_system, desc, k.first_bound, k.second_bound, ainfo, binfo, systems))

def getinfo(acc):
    Entrez.email = 'matthew.johnson@ucsf.edu'
    handle = Entrez.epost(db='protein', id=acc)
    result = Entrez.read(handle)
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
    annotations = Entrez.read(data)
    handle.close()
    return annotations

def getcolors(mylist):
    colordict = {}
    for i, item in enumerate(sum(mylist, [])):
        r = lambda: random.randint(0,255)
        color = '#%02X%02X%02X' % (r(),r(),r())
        colordict[item] = color
    return colordict

def graph(hs):
    for k,v in hs.items():
        v.regions = sorted(v.regions, key=lambda x: len(x))
        anchor = v.regions[0][0].family
        regions = [[[t.family, t.identifier] for t in x] for x in v.regions]
        colordict = getcolors( [[t.family for t in x] for x in v.regions] )
        width = len(max(regions))*20
        length = len(regions)*10
        d = draw.Drawing(width, length, origin='center', displayInline=False)
        for t, loci in enumerate(regions):
            if loci[0][0] != anchor:
                loci.reverse()
            t = t *10
            y = (length/2-10)-t
            for i, element in enumerate(loci):
                family, identifier = element
                i = i*10
                x = (-width/2)+i
                fill_opacity=0
                color = colordict[family]
                sw=0
                if identifier in island_dict:
                    sw=3
                d.append(draw.Rectangle(x,y,10,10, fill=color, stroke_width=sw, stroke_color='black'))
        d.saveSvg('{}_{}.svg'.format(k.first_bound,k.second_bound))
def cluster(regions):
    from scipy.cluster.hierarchy import fcluster, dendrogram, linkage
    #for k,v in hs.items():
    n = len(regions)
    matrix = np.zeros((n,n))
    regionsfam = [[t.family for t in x] for x in regions]
    for i, a in enumerate(regionsfam):
        for j, b in enumerate(regionsfam):
            a, b = set(a), set(b)
            u = len(a.intersection(b))/len(max([a,b]))
            matrix[i,j] = u
    plt.figure(figsize=(10, 7))
    linked = linkage(matrix, 'centroid')
    dn = dendrogram(linked)
    plt.savefig('dendo/dendo_DHS1.png')
    cut = fcluster(linked, 14, criterion="distance")
    both = list(zip(regions,cut))
    sorted_regions = sorted(both, key=lambda x: x[1])#sorted(regions, key=lambda x: len(x))  #sorted(both, key=lambda x: x[1])
    regions = sorted_regions#[x[0] for x in sorted_regions]
    print(regions)
    anchor = regions[0][0][0].family
    regions_new = [[[t.family, t.identifier] for t in x[0]] for x in regions]
    colordict = getcolors( [[t.family for t in x[0]] for x in regions] )
    width = len(max(regions_new))*20
    length = len(regions_new)*10
    d = draw.Drawing(width, length, origin='center', displayInline=False)
        #c = draw.Drawing(width, length, origin='center', displayInline=False)
    for t, loci in enumerate(regions_new):
        if loci[0][0] != anchor:
            loci.reverse()
        t = t *10
        y = (length/2-10)-t
        if True:#not any([x[1] in island_dict for x in loci]):
            for i, element in enumerate(loci):
                family, identifier = element
                i = i*10
                x = (-width/2)+i
                alpha = 1
                color = colordict[family]
        #            if identifier in island_dict:
        #                alpha = 1
        #                color = color_dict[island_dict[identifier]]
                d.append(draw.Rectangle(x,y,10,10,fill=color, fill_opacity=alpha))
                    #color = colordict[family]
                    #alpha = 1
                    #c.append(draw.Rectangle(x,y,10,10,fill=color, fill_opacity=alpha))
    d.saveSvg('hotspot_images/clustered_DHS1.svg')
        #c.saveSvg('hotspot_images/color_{}_{}.svg'.format(k.first_bound,k.second_bound))

import seaborn as sns
pal = sns.color_palette("husl", 17).as_hex()
systems = ['Kiwa', 'Wadjet', 'Hachiman', 'Thoeris', 'CBASS','Disarm','Retron','Septu', 'Shedu', 'Gabija', 'Lamassu', 'Type IV RM', 'Druantia', 'Zorya', 'Type I RM', 'Type II RM', 'Type III RM']
color_dict = dict(zip(systems,pal))
from collections import Counter
def systemcalling(all_regions):
    domains = []
    #all_regions = []
    #for k,v in hs.items():
    #    thing = ['GCA001879525qCP013113.1_APC72375.1', 'GCA001879525qCP013113.1_APC72374.1', 'GCA002194085qNIOZ01000039.1_OWI79132.1', 'GCA001879525qCP013113.1_APC72372.1', 'GCA900618275qLR130530.1_VDL06751.1']
    #    if k.first_bound in thing or k.second_bound in thing:
    #        for region in v.regions:
    #            all_regions.append(region)
    #return all_regions
    toid={}
    for region in all_regions: 
        for gene in region:
            if gene.family in toid:
                toid[gene.family].append(gene.identifier)
            else:
                toid[gene.family] = [gene.identifier,]
    tok = [[t.family for t in x] for x in all_regions]
    #remove redudnant items in list
    tokuniq ={}
    for t in tok:
        key = '|'.join(t)
        tokuniq[key] = t
    res = list(tokuniq.values())
    nets = []
    for _ in range(100):
        network = []
        [random.shuffle(x) for x in res]
        for region in res:# v.regions:
            #systempresent = False
            #if any(genefamily in island_dict for genefamily in region):
            #    systempresent = True
            if True: #systempresent == False:
                for i, gene in enumerate(region):
                    #systemname = island_dict.get(toid[gene][0], 'unk')
                    #print('{}\t{}'.format(gene,systemname))
                    if i == 0:
                        win = region[0:i+2]
                    else:
                        win = region[i:i+2]
                    if len(win) == 2:
                        x = win[0]
                        z = win[1]
                        tup = tuple(sorted([x,z]))
                        network.append(tup)
                    #d = domain_dict.get(gene, ['nf'])
                    #d = '|'.join(list(d[0]))
                    #print('{}\t{}'.format(gene,d))
                    #domains.append(d)
                    #print([domain_dict.get(x.family,'no') for x in region])
        net = dict(Counter(network))
        #for k,v in net.items():
        #    print('{}\t{}\t{}'.format(k[0],k[1],v)) 
        nets.append(net)
    stats = {}
    for net in nets:
        for k,v in net.items():
            key = '{} (interacts with) {}'.format(k[0],k[1])
            if key in stats:
                stats[key].append(v)
            else:
                stats[key] = [v,]
    for k,v in stats.items():
        mean = np.mean(v)
        sd = np.std(v)
        print('{}\t{}\t{}'.format(k,mean,sd))
    
    #print('{}\t{}\t{}'.format(k[0],k[1],v)) 
    #for k,v in calculate.items():
    #    t = total[k]
    #    print('{}\t{}'.format(k,v/t))
    #a = dict(Counter(domains))
    #b = sorted(a.items(), key = 
    #         lambda kv:(kv[1], kv[0]))
    #fig = plt.figure(figsize =(10, 7))
    #labs = [x[0] for x in b]
    #data = [x[1] for x in b]
    #plt.pie(data, labels = labs, autopct='%.2f%%',pctdistance=0.9)
    #plt.savefig('pie.png')
if __name__ == "__main__":
    print('getting rep dict')
    cluster_dict = get_cluster_dict() #returns a dictionary of the mmseqs cluster representatives
    print('getting contig protein dict')
    contig_protein_dict, db_size = get_contig_dict(cluster_dict)
    print(db_size)
    #db_size = 3239
    #pickle.dump(cluster_dict, open( "DHS_cluster_dict.p", "wb" ) )
    #pickle.dump(contig_protein_dict, open( "DHS_contig_dict.p", "wb" ) )
    #cluster_dict = pickle.load( open( "DHS_cluster_dict.p", "rb" ) )
    #contig_protein_dict = pickle.load( open( "DHS_contig_dict.p", "rb" ) )
    print('running get_regions2()')
    regions_dict = get_regions2(contig_protein_dict)
    #pickle.dump(regions_dict, open( "regions_dict_062021.p", "wb" ))
    #print('loading regions_dict')
    #regions_dict = pickle.load( open( "regions_dict_062021.p", "rb" ) )
    print('done') 
    #print('reading island data')
    island_dict = read_ISLAND()
    #pickle.dump(island_dict, open( "island_dict_062221.p", "wb" ) )
    #print('done')
    print('running DHS()')
    spots = DHS(regions_dict, island_dict, db_size, 0.8)
    i+1
    #pickle.dump(spots, open( "spots_062021_2.p", "wb" ) )
    #island_dict = pickle.load( open("island_dict_062221.p", 'rb' ) )
    #spots = pickle.load( open("spots_062021_2.p", 'rb' ) )
    for k,v in spots.items():
        allsystems = set()
        DHS = dict()
        for locus in v.regions:
            for gene in locus:
                id  = gene.family.split('_')[1]
                if id in island_dict:
                    if island_dict[id] in DHS:
                        DHS[island_dict[id]].append(locus)
                    else:
                        DHS[island_dict[id]] = [locus, ]
                try:
                    allsystems.add(island_dict[ id ])
                except:
                    continue
        if len(allsystems) > 9:
            print([k.first_bound, k.second_bound])
            for k,v in DHS.items():
                print(k)
                for locus in v:
                    print(locus[0].contig)
                    for gene in locus:
                        id  = gene.family.split('_')[1]
                        if id in island_dict:
                            system = island_dict[id]
                            print([system, gene.identifier])

    ##all_DHS1 = pickle.load( open('all_DHS1.p','rb'))

    ##domain_dict = pickle.load( open('candidte_DHS_unq/DHS_domain_regions.p', 'rb') )
    ##for x in all_DHS1:
        #print(x[0].contig)
    ##    for y in x:
    ##        print([y.family, domain_dict.get(y.family,'nf'), '|', island_dict.get(y.identifier, 'ns')]) 
    #DHS2(hs)
    #cluster(all_DHS1)
    #systemcalling(all_DHS1)
    #pickle.dump(all_DHS1, open( "all_DHS1.p", "wb" ) )
## incorporate ISLAND annotation into SeqUnit
## count how many times boundary pairs sandwich a DP
## di


