import os
import subprocess
import pandas as pd
import pysam
from collections import namedtuple, OrderedDict
from sets import Set
from progressbar import AnimatedMarker, Bar, BouncingBar, Counter, ETA, \
    ProgressBar, ReverseBar, RotatingMarker, \
    SimpleProgress, Timer


def generate_index(bam_input_filename, index_filename, three_prime):

    
        
    number_of_reads = int(subprocess.check_output(["samtools", "view", "-c", bam_input_filename]))

    widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
    pbar = ProgressBar(widgets=widgets, maxval=number_of_reads).start()
    i = 0  

    GeneInfo = namedtuple('GeneInfo', ['start_pos', 'end_pos'])

    gene_dic = OrderedDict([('rrnA', GeneInfo(4035153, 4040906)), 
                            ('rrnB', GeneInfo(4165428, 4172057)),
                            ('rrnC', GeneInfo(3941327, 3946872)),
                            ('rrnD', GeneInfo(3423194, 3429236)),
                            ('rrnE', GeneInfo(4207532, 4213234)),
                            ('rrnG', GeneInfo(2725746, 2731600)),
                            ('rrnH', GeneInfo(223408, 229167))])

    columns = gene_dic.keys()

    samfile = pysam.Samfile(bam_input_filename, 'rb')
    names = set()
    for read in samfile:
        names.add(read.qname)
        i += 1
        pbar.update(i)
    pbar.finish()
         
    names_list = list(names)

    df = pd.DataFrame(index=names_list, columns=columns)

    samfile = pysam.Samfile(bam_input_filename, 'rb')
               
    widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
    pbar = ProgressBar(widgets=widgets, maxval=number_of_reads).start()
    i = 0

    
    for read in samfile:
        
        if three_prime == True:
            forward_pos = read.pos + read.rlen
            reverse_pos = read.pos + 1
        else: 
            forward_pos = read.pos + 1
            reverse_pos = read.pos + read.rlen
            
        i += 1
        pbar.update(i)
        if read.is_unmapped == False:
            if read.is_reverse == False:
                position = forward_pos
                position_marker = forward_pos
            elif read.is_reverse == True:
                position = reverse_pos
                position_marker = '-' + str(reverse_pos)
            for gene, gene_data in gene_dic.items():
                if (position > gene_data.start_pos and position < gene_data.end_pos):
                    df[gene][read.qname] = position_marker

                        
    pbar.finish()
    samfile.close()
    df.index.name = 'Readname'
    df.reset_index(level=0, inplace=True)
    df.to_csv(index_filename, index=True ,header=True)

    '''index_filename_2 = index_filename.replace('_index.csv', '_readname.csv' )
    del df['Readname']
    df.to_csv(index_filename_2, index=True, header= true)'''
    

