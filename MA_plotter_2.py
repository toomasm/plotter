import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import pandas as pd
from collections import namedtuple, OrderedDict, Counter, defaultdict
import math
import numpy as np
from Bio import SeqIO

import hdf5_gen



def by_strand_sorter(index, add_to_pos_reads, add_to_neg_reads, counter_pos, counter_neg):
    if index in counter_pos and index in counter_neg:
        add_to_pos_reads.append(counter_pos[index])
        add_to_neg_reads.append(counter_neg[index])
    elif index not in counter_pos and index in counter_neg:
        add_to_neg_reads.append(counter_neg[index])
        add_to_pos_reads.append(1)
    elif index in counter_pos and index not in counter_neg:
        add_to_pos_reads.append(counter_pos[index])
        add_to_neg_reads.append(1)
    else:
        add_to_pos_reads.append(1)
        add_to_neg_reads.append(1)

        

def data_generator(index, selected_column, standard_colour, ACA_colour, MqsR_colour, marker_symbol):

    GeneInfo = namedtuple('GeneInfo', ['start_pos', 'end_pos', 'column_label'])
    plot_data = namedtuple('PlottingData', ['nucl_data' ,'y_pos', 'y_neg' ,'colour', 'symbol'])

    gene_dic = OrderedDict([('rrnA', GeneInfo(4035153, 4040906, 'rrnA')),
                            ('rrnB', GeneInfo(4165428, 4172057, 'rrnB')),
                            ('rrnC', GeneInfo(3941327, 3946872, 'rrnC')),
                            ('rrnD', GeneInfo(3423194, 3429236, 'rrnD')),
                            ('rrnE', GeneInfo(4207532, 4213234, 'rrnE')),
                            ('rrnG', GeneInfo(2725746, 2731600, 'rrnG')),
                            ('rrnH', GeneInfo(223408, 229167, 'rrnH'))])

    all_operons = ['rrnA', 'rrnB', 'rrnC', 'rrnD', 'rrnE', 'rrnG', 'rrnH']
    operon_positions = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', 'P7']

    ref_genome_fasta = 'genbank_mg1655.fasta'

    #selected_column = 'rrnA'
    selected_positions_column = gene_dic[selected_column].column_label

    
    dataframe = pd.read_csv(index)
    values_list_pos = []
    values_list_neg = []
    for value in  dataframe[selected_positions_column].dropna():
        if '+AC0-' in str(value):                
            value = str(value)[4:]
            values_list_pos.append(int(float(value)))
        elif '-' in str(value):
            value = str(value)[1:]
            values_list_neg.append(int(float(value)))
        else:
            values_list_pos.append(int(float(value)))

    #deletes dataframe to free memory
    del dataframe
    
    #Counts how many reads are mapped to every position on both pos and neg strand.
    #Returns dict
    counter_pos = Counter(values_list_pos)
    counter_neg = Counter(values_list_neg)
   
    #Create empty lists for storing reads 5 prime position count on genome.
    # nucl_data = genome position, y = read count, pos and neg determine the strand orientation
    #colour = separates cutting seq from rest by colour, stores a color of dot on a plot for each nucleotide position
    #in appropriate colour. No cutting site = blue.
    
    nucl_data = []
    y_pos = []
    y_neg = []
    colour = []
    
    for fasta in SeqIO.parse(ref_genome_fasta, "fasta"):
        pass
    
    for index in range(gene_dic[selected_column].start_pos, gene_dic[selected_column].end_pos):

        if str(fasta.seq[index-1: index+2]) == 'ACA':
            colour.append(ACA_colour)
        elif str(fasta.seq[index-2: index+1]) == 'GCT' or str(fasta.seq[index-2: index+1]) == 'GCA':
            colour.append(MqsR_colour)
        else:
            colour.append(standard_colour)
        nucl_data.append(index)
        by_strand_sorter(index, y_pos, y_neg, counter_pos, counter_neg)

    symbol = marker_symbol    
                           
    return plot_data(nucl_data, y_pos, y_neg, colour, symbol)

def make_plot(_input_, selected_column, scatter_flag, MA_flag, three_prime_flag):
    """ Creates a MA plot or scatter plot against chosen rrnA operon from two CSV index files (index, positions, mapped against how many).

        Scatter plot can be made from up to 4 files. 3 files can be plotted againist each other normaly. With 4 files the program assumes
        that you are looking at 3' seq reads that have been polyA trimmed. Two first ones are treated as ok and shady reads of the same sample
        and treatment and the 3-rd and fourth as ok and shady reads of second sample.

        Because my system has memory limitations with pandas I strip down the indexes before plotting them here. All reads data
        I want to show is present in the CSV file."""

    #Checks whether three_prime seq flag has been set. If yes makes symbols for ok reads 'o' and for shady reads '' 
    hdf_filename = hdf5_gen.make_hdf_filename(_input_, scatter_flag, MA_flag, three_prime_flag, selected_column)

    try:
        data_dic = hdf5_gen.get_hdf_data(hdf_filename)

    except IOError:

        if three_prime_flag:
            symbol_1 = 'o'
            symbol_2 = 'D'

        else:
            symbol_1 = 'o'
            symbol_2 = 'o'
        
        #looks how many input files have been given
        input_number = len(_input_)

        data_dic = {}
        data_dic['data1'] = data_generator(_input_[0], selected_column, 'b', 'r', 'm', symbol_1)

        #Creates data necessary for plotting for every inputfile.
        if input_number >= 2:
            data_dic['data2'] = data_generator(_input_[1], selected_column, 'b', 'r', 'm', symbol_2)
            
        if input_number >= 3:
            data_dic['data3'] = data_generator(_input_[2], selected_column, 'g', 'c', 'm', symbol_1)
            
        if input_number == 4:
            data_dic['data4'] = data_generator(_input_[3], selected_column, 'g', 'c', 'm', symbol_2)
        
        hdf5_gen.add_to_hdf(data_dic, hdf_filename)

    def onpick(event):
        
        #Gets the index of datapoint (index of X and Y)
        ind = event.ind

        #Retrieves information from lists based on index of the datapoint.
        #The information to be displayed for datapoint has the same index as datapint.

        nucleotide_pos = np.take(data_dic['data1'].nucl_data, ind)
        sample_1_pos_read_count = np.take(data_dic['data1'].y_pos, ind)
        sample_2_pos_read_count = np.take(data_dic['data2'].y_pos, ind)
        sample_1_neg_read_count = np.take(data_dic['data1'].y_neg, ind)
        sample_2_neg_read_count = np.take(data_dic['data2'].y_neg, ind)

        if input_number >= 3:
            sample_3_pos_read_count = np.take(data_dic['data3'].y_pos, ind)
            sample_3_neg_read_count = np.take(data_dic['data3'].y_neg, ind)
        if input_number == 4:
            sample_4_pos_read_count = np.take(data_dic['data4'].y_pos, ind)
            sample_4_neg_read_count = np.take(data_dic['data4'].y_neg, ind)

        #Prints out data for datapoints in nucleotide position (X-axis)  that was covered with a click.
        #Prints out the read counts for + and - strand for nucleotide position in each sample present.
        for array_ind in range(len(nucleotide_pos)):
           
            if input_number == 3:    
                print ("Nucleotide position: {0} , s1 pos: {1} , s2 pos: {2}, s3 pos: {3} \n \
                             s1 neg: {4} , s2 neg: {5}, s3 neg: {6}".format(nucleotide_pos[array_ind], sample_1_pos_read_count[array_ind], 
                                                                            sample_2_pos_read_count[array_ind], sample_3_pos_read_count[array_ind],
                                                                            sample_1_neg_read_count[array_ind],sample_2_neg_read_count[array_ind],
                                                                            sample_3_neg_read_count[array_ind]))

            elif input_number == 4:
                print ("Nucleotide position: {0} , s1 pos: {1} , s2 pos: {2}, s3 pos: {3}, s4 pos: {4} \n \
                              s1 neg: {5} , s2 neg: {6}, s3 neg: {7}, s4 neg: {8}".format(nucleotide_pos[array_ind], sample_1_pos_read_count[array_ind], 
                                                                               sample_2_pos_read_count[array_ind], sample_3_pos_read_count[array_ind],
                                                                               sample_4_pos_read_count[array_ind], sample_1_neg_read_count[array_ind],
                                                                               sample_2_neg_read_count[array_ind], sample_3_neg_read_count[array_ind],
                                                                               sample_4_neg_read_count[array_ind]))

            else:
                
                print ("Nucleotide position: {0} , s1 pos: {1} , s2 pos {2} \n \
                            s1 neg: {3} , s2 neg {4}".format(nucleotide_pos[array_ind], sample_1_pos_read_count[array_ind], 
                                                             sample_2_pos_read_count[array_ind], sample_1_neg_read_count[array_ind],
                                                             sample_2_neg_read_count[array_ind]))
                
    #Creates a MA plot. X is A and Y is M (look a bit back and you see they represent).
    if MA_flag:

        #Generates data for MA plot.
        MA_X = [(math.log(float(y1), 2) + math.log(float(y2), 2))/2 for y1, y2 in zip(data_dic['data1'].y_pos, data_dic['data2'].y_pos)]
        MA_Y = [math.log((float(y1)/y2), 2) for y1, y2 in zip(data_dic['data1'].y_pos, data_dic['data2'].y_pos)]
        
        fig1, ax = plt.subplots()
        ax.set_ylabel('M', fontsize=25)
        ax.set_xlabel('A', fontsize=25)
        ax.scatter(MA_X, MA_Y, alpha=0.5, c=data_dic['data1'].colour, linewidths=( 0, 0, 0), picker=True)
        fig1.canvas.mpl_connect('pick_event', onpick)

        fig1.savefig('test5.pdf')
        
    
    #Creates a scatter plot. X is nucleotide position (for selected rrnA operon), Y is read count on specific position.
    if scatter_flag:
        
        fig2, ax = plt.subplots()
        ax.set_ylabel('Read Counts', fontsize=25)
        ax.set_xlabel('Nucleotide Position', fontsize=25)

        for data_key, data_nt in data_dic.items():#[data1, data2, data3, data4]:
            ax.scatter(data_nt.nucl_data, data_nt.y_pos, alpha=0.5, c=data_nt.colour, linewidths=( 0, 0, 0), picker=True, marker = data_nt.symbol)
            ax.scatter(data_nt.nucl_data, [-1 * data for data in data_nt.y_neg], alpha=0.5, c=data_nt.colour, linewidths=( 0, 0, 0), picker=True, marker = data_nt.symbol)
            
        if three_prime_flag:
            fig2.canvas.mpl_connect('pick_event', onpick)
        else:
            fig2.canvas.mpl_connect('pick_event', onpick)
            
        fig2.savefig('test6.pdf')

    plt.show()