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
    """ Generates index file that has read names column and rRNA operons columns. If a read is mapped to some rRNA operon the mapped
    position in that operon will be added to the right operon column for that read. One read can be mapped to several operons. E. coli
    has 7 rRNA operons.
    
    Output is csv file readable by pandas."""
    
        
    number_of_reads = int(subprocess.check_output(["samtools", "view", "-c", bam_input_filename]))

    widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
    pbar = ProgressBar(widgets=widgets, maxval=number_of_reads).start()

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
    
    #Make a set of read names. The issue is that we have same reads mapped to different locations and are represented as
    #different reads in bam file. So for constructing the index file we will need one read name only once.
    names = set()
    for index, read in enumerate(samfile):
        names.add(read.qname)
        pbar.update(index)
    pbar.finish()
         
    names_list = list(names)

    #Generate dictionary for read data. 
    data_dic = {}

    #Dict keys shall be read names.
    for n in names:
        data_dic[n] = ["","","","","","",""]

    #Delete set containing read names.
    del names
    print len(data_dic)
    
    samfile = pysam.Samfile(bam_input_filename, 'rb')
               
    widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]

    print('number_of_reads', number_of_reads)
    pbar = ProgressBar(widgets=widgets, maxval=number_of_reads).start()

    for index, read in enumerate(samfile):
        pbar.update(index)
                
        if read.is_unmapped == False:
            if three_prime == True:
                forward_pos = read.pos + read.rlen
                reverse_pos = read.pos + 1
            else: 
                forward_pos = read.pos + 1
                reverse_pos = read.pos + read.rlen

            if not read.is_reverse:
                position = forward_pos
                position_marker = forward_pos
            else:
                position = reverse_pos
                position_marker = '-' + str(reverse_pos)
            n= 0
            for gene, gene_data in gene_dic.items():
                if (position > gene_data.start_pos and position < gene_data.end_pos):
                    data_dic[read.qname][n] = position_marker
                n += 1

    df = pd.DataFrame.from_dict(data_dic, orient='index', dtype=None)
    df.rename(columns={0: columns[0], 1: columns[1], 2: columns[2], 3: columns[3],4: columns[4],
                       5: columns[5], 6: columns[6]}, inplace=True)
                    
    pbar.finish()
    samfile.close()
    df.index.name = 'Readname'
    df.reset_index(level=0, inplace=True)
    df.to_csv(index_filename, index=True ,header=True)

    '''index_filename_2 = index_filename.replace('_index.csv', '_readname.csv' )
    del df['Readname']
    df.to_csv(index_filename_2, index=True, header= true)'''
    

