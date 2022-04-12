#pipeline
#build-icm genome.icm < input_fna
#glimmer3 input_fna genome.icm genome
#extract input_fna genome.predict > genome.fna
from Bio import SeqIO
from Bio.Seq import Seq
import os
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import subprocess
import random
import string
import glob
import sys
def randomString(stringLength=8):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

def makedir(working_directory):
    f  =  os.path.join(working_directory, randomString(6))
    os.mkdir(f)
    directories = ['tmp']
    for basename in directories:
        os.mkdir( os.path.join(f, basename) )
    return f

def build_icm(filename, tmp_dir):
    print('building ICM')
    #filename = os.path.join(tmp_dir, 'user_fasta.fna')
    output_dir = os.path.join(tmp_dir, 'genome.icm')
    subprocess.run(['/bin/build-icm ' + output_dir + ' < ' + filename], shell=True)    

def glimmer3(filename, tmp_dir):
    print('running glimmer3')
    #filename_fna = os.path.join(tmp_dir, 'user_fasta.fna')
    genome_icm = os.path.join(tmp_dir, 'genome.icm')
    subprocess.run(['/bin/glimmer3 -l --gc_percent 67.2 ' + filename + ' ' + genome_icm + ' ' + os.path.join(tmp_dir,'tmp','genome')], shell=True)

def extract(filename, tmp_dir):
    print('extracting')
    #filename_fna = os.path.join(tmp_dir, 'user_fasta.fna')
    subprocess.run(['/bin/extract ' + filename + ' ' + os.path.join(tmp_dir, 'tmp', 'genome.predict') + ' > ' + os.path.join(tmp_dir, 'tmp', 'genome.fna')], shell=True)

def makeGenbank(filename_fna, tmp_dir, information):
    print('making genbank')
    #fn = os.path.basename(os.path.normpath(filename_fna))
    #name = fn.split('.')[0]
    name = 'output.gbk'
    #filename_fna = os.path.join(tmp_dir, 'user_fasta.fna')
    with open(filename_fna) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            fna = record.seq
            print(len(fna))
    #create seqRecord
    record_to_write = SeqRecord(Seq(str(fna)), id="Personal", description='Personal', name=name, annotations=information)
    #create features
    filename_predict = os.path.join(tmp_dir, 'tmp', 'genome.fna')
    with open(filename_predict) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            print(record.seq)
            print(record.seq.translate())
            qualifiers = {
            'locus_tag': record.id,
            'product': 'unknown',
            'protein_id': record.id,
            'translation': record.seq.translate(table=11)
            }
            
            start, end = record.description.split()[1:3]
            start = int(start)
            end = int(end)
            print(record.description)
            print([start,end])
            strand = 1
            if end < start:
                start, end = end, start
                strand = -1
            new_one = SeqFeature(location = FeatureLocation(start=int(start),
                                                 end=int(end),
                                                 strand=strand),
                                                 type='CDS',
                                                 qualifiers=qualifiers)
            record_to_write.features.append(new_one)
    #fn = os.path.basename(os.path.normpath(filename_fna))
    #output_filename = fn.split('.')[0] + '.gbk'
    handle = open( os.path.join(tmp_dir, 'tmp', 'output.gbk'), 'w')
    SeqIO.write(record_to_write, handle, 'genbank')
    handle.close()
    #return os.path.join(tmp_dir, 'tmp', output_filename)

def fasta2Gembase(files):
    information = {
        'molecule_type':'DNA',
        'topology': 'linear',
        'data_file_division': 'BCT',
        'date': '2-JUNE-2021',
        'accessions': ['PERSONAL'],
        'sequence_version': 1,
        'source': 'Phage',
        'organism': 'Pseudomonas phage',
        'taxonomy': ['Bacteria',
          'Proteobacteria',
          'Gammaproteobacteria',
          'Pseudomonadales',
          'Pseudomonadaceae',
          'Pseudomonas'],
        }
    #for filename_fna in files:
    filename_fna = sys.argv[1]
    tmp_dir = makedir('/hdd-roo/DHS/meta_gbk/')
    build_icm(filename_fna, tmp_dir)
    glimmer3(filename_fna, tmp_dir)
    extract(filename_fna, tmp_dir)
    output_filename = makeGenbank(filename_fna, tmp_dir, information)
    print('output at {}'.format(output_filename))

fasta2Gembase(glob.glob('/hdd-roo/DHS/meta_fna/*fna'))
