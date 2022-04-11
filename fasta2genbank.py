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

def randomString(stringLength=8):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

def makedir(working_directory):
    f  =  os.path.join(working_directory, randomString(4))
    os.mkdir(f)
    directories = ['tmp']
    for basename in directories:
        os.mkdir( os.path.join(f, basename) )
    return f

def build_icm(tmp_dir):
    filename = os.path.join(tmp_dir, 'user_fasta.fna')
    output_dir = os.path.join(tmp_dir, 'genome.icm')
    subprocess.run(['/bin/build-icm ' + output_dir + ' < ' + filename], shell=True)    

def glimmer3(tmp_dir):
    filename_fna = os.path.join(tmp_dir, 'user_fasta.fna')
    genome_icm = os.path.join(tmp_dir, 'genome.icm')
    subprocess.run(['/bin/glimmer3 ' + filename_fna + ' ' + genome_icm + ' ' + os.path.join(tmp_dir,'tmp','genome')], shell=True)

def extract(tmp_dir):
    filename_fna = os.path.join(tmp_dir, 'user_fasta.fna')
    subprocess.run(['/bin/extract ' + filename_fna + ' ' + os.path.join(tmp_dir, 'tmp', 'genome.predict') + ' > ' + os.path.join(tmp_dir, 'tmp', 'genome.fna')], shell=True)

def makeGenbank(tmp_dir, information):
    #fn = os.path.basename(os.path.normpath(filename_fna))
    #name = fn.split('.')[0]
    name = 'output.gbk'
    filename_fna = os.path.join(tmp_dir, 'user_fasta.fna')
    with open(filename_fna) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            fna = record.seq
            print(len(fna))
    #create seqRecord
    record_to_write = SeqRecord(Seq(str(fna), generic_dna), id="Personal", description='Personal', name=name, annotations=information)
    #create features
    filename_predict = os.path.join(tmp_dir, 'tmp', 'genome.fna')
    with open(filename_predict) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            qualifiers = {
            'locus_tag': record.id,
            'product': 'unknown',
            'protein_id': record.id,
            'translation': record.seq.translate(table=11)
            }
            
            start, end = record.description.split()[1:3]
            start = int(start)
            end = int(end)
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

if __name__ == '__main__':
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

    filename_fna = '/home/jbd_apps/fasta_to_genbank/KMV_phages/gdFB465.fasta'
    tmp_dir = makedir('/home/jbd_apps/fasta_to_genbank')
    build_icm(filename_fna, tmp_dir)
    glimmer3(filename_fna, tmp_dir)
    extract(filename_fna, tmp_dir)
    output_filename = makeGenbank(tmp_dir, information, filename_fna)
    print('output at {}'.format(output_filename)) 
