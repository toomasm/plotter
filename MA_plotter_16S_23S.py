import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import pandas as pd
from collections import namedtuple, OrderedDict, Counter, defaultdict
import math
import numpy as np
from Bio import SeqIO
import hdf5_gen
from progressbar import AnimatedMarker, Bar, BouncingBar, ETA, \
    ProgressBar, ReverseBar, \
    SimpleProgress, Timer

LENGTH_OF_MATURE_16S = 1542
#LENGTH_OF_MATURE_23S = 2905
LENGTH_OF_MATURE_23S = 2925

def mature_ribosome_RNA_divider(counter, rRNA):
    """ Takes dictionary containing all rRNA operon positions as keys and read counts in position as values as input.
    Generates a new dictionary representeing a blank 16S or 23S rRNA. Adds all those reads in appropriate position on
    the blank rRNA."""
    
    ribosome_dict = {}
    if rRNA == "16S":
        subunit = 0
    elif rRNA == "23S":
        subunit = 1
        
    for operon in all_operons:
        for value in range (gene_dic[operon][subunit].start_pos, gene_dic[operon][subunit].end_pos + 1):
            if value in counter:
                new_value = value - gene_dic[operon][subunit].start_pos + 1
                if new_value not in ribosome_dict:
                    ribosome_dict[new_value] = counter[value]
                else:
                    ribosome_dict[new_value] = ribosome_dict[new_value] + counter[value] + 1
                    
        
    return ribosome_dict


    
def by_strand_sorter(index, add_to_pos_reads, add_to_neg_reads, counter_pos, counter_neg):
    """ Adds read count information to positive and negative strand lists. If read count
    in position = 0 add 1 instead (MA plot cant be made with 0 values)."""
    
    if index in counter_pos:
        add_to_pos_reads.append(counter_pos[index])
    else:
        add_to_pos_reads.append(1)
        
    if index in counter_neg:
        add_to_neg_reads.append(counter_neg[index])
    else:
        add_to_neg_reads.append(1)

plot_data = namedtuple('PlottingData', ['nucl_data_16S', 'nucl_data_23S' ,'y_pos_16S', 'y_neg_16S', 'y_pos_23S',
                                        'y_neg_23S', 'colour_16S', 'colour_23S', 'symbol'])
        
def data_generator(index, standard_colour, ACA_colour, MqsR_colour, symbol):

    global all_operons, gene_dic
    
    GeneInfo = namedtuple('GeneInfo', ['start_pos', 'end_pos', 'column_label'])

    
    _16S = namedtuple('_16S', ['start_pos', 'end_pos'])
    _23S = namedtuple('_23S', ['start_pos', 'end_pos'])

    '''gene_dic = OrderedDict([('rrnA', (_16S(4035531, 4037072), _23S(4037519, 4040423))),
                            ('rrnB', (_16S(4166659, 4168200), _23S(4168641, 4171544))),
                            ('rrnC', (_16S(3941808, 3943349), _23S(3943704, 3946607))),
                            ('rrnD', (_16S(3427221, 3428762), _23S(3423880, 3426783))),
                            ('rrnE', (_16S(4208147, 4209688), _23S(4210043, 4212946))),
                            ('rrnG', (_16S(2729616, 2731157), _23S(2726281, 2729184))),
                            ('rrnH', (_16S(223771, 225312), _23S(225759, 228662)))])'''

    gene_dic = OrderedDict([('rrnA', (_16S(4035531, 4037072), _23S(4037499, 4040423))),
                            ('rrnB', (_16S(4166659, 4168200), _23S(4168621, 4171544))),
                            ('rrnC', (_16S(3941808, 3943349), _23S(3943684, 3946607))),
                            ('rrnD', (_16S(3427221, 3428762), _23S(3423880, 3426803))),
                            ('rrnE', (_16S(4208147, 4209688), _23S(4210023, 4212946))),
                            ('rrnG', (_16S(2729616, 2731157), _23S(2726281, 2729204))),
                            ('rrnH', (_16S(223771, 225312), _23S(225759, 228662)))])

        
    all_operons = ['rrnA', 'rrnB', 'rrnC', 'rrnD', 'rrnE', 'rrnG', 'rrnH']
    
    ref_genome_fasta_16S = 'rrsA.fasta'
    #ref_genome_fasta_23S = 'rrlA.fasta'
    ref_genome_fasta_23S = 'rrlA_extended.fasta'
    
    dataframe = pd.read_csv(index)

    values_list_pos = []
    values_list_neg = []


    
    for operon in all_operons:
        
        number_of_reads = len(dataframe)
        widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
        pbar = ProgressBar(widgets=widgets, maxval=number_of_reads).start()
        i = 0

        for value in dataframe[operon].dropna():

            #if value in range (1, 1453):
            if '+AC0-' in str(value):
                value = int(float(str(value)[4:]))
                #if gene_dic[selected_column][0].start_pos <= value <= gene_dic[selected_column][0].end_pos:
                values_list_pos.append(value)
            if '-' in str(value):
                value = int(float(str(value)[1:]))
                #if gene_dic[selected_column][0].start_pos <= value <= gene_dic[selected_column][0].end_pos:
                values_list_neg.append(value)             
            else:
                #if value in range (gene_dic[selected_column][0].start_pos, gene_dic[selected_column][0].end_pos):
                values_list_pos.append(int(float(value)))           
            
            
            i += 1
            pbar.update(i)
       
        pbar.finish()
        print len(dataframe)
        dataframe_2 = dataframe[dataframe[operon].isnull()]
        #dataframe = dataframe.dropna(subset=[operon]) 
        dataframe = dataframe_2
        print len(dataframe)

        
    #deletes dataframes to free memory
    del dataframe
    del dataframe_2

    # Creates a dictionary where key is position and value is count of reads at that position.
    counter_pos = Counter(values_list_pos)
    counter_neg = Counter(values_list_neg)

    # Puts all the reads from different operons on the same blank operon for both 16S and 23S mature RNA-s.
    counter_pos_16S = mature_ribosome_RNA_divider(counter_pos, "16S")
    counter_pos_23S = mature_ribosome_RNA_divider(counter_pos, "23S")
    counter_neg_16S = mature_ribosome_RNA_divider(counter_neg, "16S")
    counter_neg_23S = mature_ribosome_RNA_divider(counter_neg, "23S")
        
        
    #Create empty lists for storing reads 5 pr 3 prime position count on genome.
    #nucl_data = genome position, y = read count, pos and neg determine the strand orientation
    #colour = separates cutting seq from rest by colour, stores a color of dot on a plot for each nucleotide position
    #in appropriate colour. No cutting site = blue.
    
    nucl_data_16S = []
    nucl_data_23S = []
    y_pos_16S = []
    y_neg_16S = []
    y_pos_23S = []
    y_neg_23S = []
    colour_16S = []
    colour_23S = []
    
    for fasta in SeqIO.parse(ref_genome_fasta_16S, "fasta"):
        pass
    
    for index in range(1, LENGTH_OF_MATURE_16S + 1):

        if str(fasta.seq[index-1: index+2]) == 'ACA':
            colour_16S.append(ACA_colour)
        elif str(fasta.seq[index-2: index+1]) == 'GCT' or str(fasta.seq[index-2: index+1]) == 'GCA':
            colour_16S.append(MqsR_colour)
        else:
            colour_16S.append(standard_colour)
        nucl_data_16S.append(index)
        by_strand_sorter(index, y_pos_16S, y_neg_16S, counter_pos_16S, counter_neg_16S)

    for fasta in SeqIO.parse(ref_genome_fasta_23S, "fasta"):
        pass

    for index in range(1, LENGTH_OF_MATURE_23S + 1):

        if str(fasta.seq[index-1: index+2]) == 'ACA':
            colour_23S.append(ACA_colour)
        elif str(fasta.seq[index-2: index+1]) == 'GCT' or  str(fasta.seq[index-2: index+1]) == 'GCA':
            colour_23S.append(MqsR_colour)  
        else:
            colour_23S.append(standard_colour)
        nucl_data_23S.append(index-20)
        by_strand_sorter(index, y_pos_23S, y_neg_23S, counter_pos_23S, counter_neg_23S)

    
    return plot_data(nucl_data_16S, nucl_data_23S, y_pos_16S, y_neg_16S, y_pos_23S, y_neg_23S, colour_16S, colour_23S, symbol)
   
                           

def make_plot_all_reads(_input_, scatter_flag, MA_flag, three_prime_flag):
    """ Creates a MA plot or scatter plot against chosen rrnA operon from two CSV index files (index, positions, mapped against how many).
        Because my system has memory limitations with pandas I strip down the indexes before plotting them here. All reads data
        I want to show is present in the CSV file."""

    hdf_filename = hdf5_gen.make_hdf_filename(_input_, scatter_flag, MA_flag, three_prime_flag)

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
        data_dic['data1'] = data_generator(_input_[0], 'b', 'r', 'm', symbol_1)
        #Creates data necessary for plotting for every inputfile.

        if input_number >= 2:
            data_dic['data2'] = data_generator(_input_[1], 'g', 'r', 'm', symbol_2)
            
        if input_number >= 3:
            data_dic['data3'] = data_generator(_input_[2], 'y', 'c', 'm', symbol_1)
            
        if input_number == 4:
            data_dic['data4'] = data_generator(_input_[3], 'g', 'c', 'm', symbol_2)

        hdf5_gen.add_to_hdf(data_dic, hdf_filename)
 
    def onpick_16S(event):
        
        #Gets the index of datapoint (index of X and Y)
        ind = event.ind

        #Retrieves information from lists based on index of the datapoint.
        #The information to be displayed for datapoint has the same index as datapint.

        nucleotide_pos = np.take(data_dic['data1'].nucl_data_16S, ind)
        sample_1_pos_read_count = np.take(data_dic['data1'].y_pos_16S, ind)
        sample_2_pos_read_count = np.take(data_dic['data2'].y_pos_16S, ind)
        sample_1_neg_read_count = np.take(data_dic['data1'].y_neg_16S, ind)
        sample_2_neg_read_count = np.take(data_dic['data2'].y_neg_16S, ind)

        if input_number >= 3:
            sample_3_pos_read_count = np.take(data_dic['data3'].y_pos_16S, ind)
            sample_3_neg_read_count = np.take(data_dic['data3'].y_neg_16S, ind)
        if input_number == 4:
            sample_4_pos_read_count = np.take(data_dic['data4'].y_pos_16S, ind)
            sample_4_neg_read_count = np.take(data_dic['data4'].y_neg_16S, ind)

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
                

    def onpick_23S(event):
        
        #Gets the index of datapoint (index of X and Y)
        ind = event.ind

        #Retrieves information from lists based on index of the datapoint.
        #The information to be displayed for datapoint has the same index as datapint.

        nucleotide_pos = np.take(data_dic['data1'].nucl_data_23S, ind)
        sample_1_pos_read_count = np.take(data_dic['data1'].y_pos_23S, ind)
        sample_2_pos_read_count = np.take(data_dic['data2'].y_pos_23S, ind)
        sample_1_neg_read_count = np.take(data_dic['data1'].y_neg_23S, ind)
        sample_2_neg_read_count = np.take(data_dic['data2'].y_neg_23S, ind)

        if input_number >= 3:
            sample_3_pos_read_count = np.take(data_dic['data3'].y_pos_23S, ind)
            sample_3_neg_read_count = np.take(data_dic['data3'].y_neg_23S, ind)
        if input_number == 4:
            sample_4_pos_read_count = np.take(data_dic['data4'].y_pos_23S, ind)
            sample_4_neg_read_count = np.take(data_dic['data4'].y_neg_23S, ind)

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

    if MA_flag and scatter_flag:
        fig1, (ax11, ax12) = plt.subplots(2)
        fig2, (ax21, ax22) = plt.subplots(2)            
    else:
        fig1, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots()
    
    if MA_flag:
        if scatter_flag:
            ax1 = ax11
            ax2 = ax21
            
        #Generates data for MA plot.
        MA_X_16S = [(math.log(float(y1), 2) + math.log(float(y2), 2))/2 for y1, y2 in zip(data_dic['data1'].y_pos_16S, data_dic['data2'].y_pos_16S)]
        MA_Y_16S = [math.log((float(y1)/y2), 2) for y1, y2 in zip(data_dic['data1'].y_pos_16S, data_dic['data2'].y_pos_16S)]

        MA_X_23S = [(math.log(float(y1), 2) + math.log(float(y2), 2))/2 for y1, y2 in zip(data_dic['data1'].y_pos_23S, data_dic['data2'].y_pos_23S)]
        MA_Y_23S = [math.log((float(y1)/y2), 2) for y1, y2 in zip(data_dic['data1'].y_pos_23S, data_dic['data2'].y_pos_23S)]

        
        ax1.set_ylabel('M', fontsize=20)
        ax1.set_xlabel('A', fontsize=20)
        ax1.scatter(MA_X_16S, MA_Y_16S, alpha=0.5, c=data_dic["data1"].colour_16S, linewidths=( 0, 0, 0), picker=True)

        ax2.set_ylabel('M', fontsize=20)
        ax2.set_xlabel('A', fontsize=20)
        ax2.scatter(MA_X_23S, MA_Y_23S, alpha=0.5, c=data_dic["data1"].colour_23S, linewidths=( 0, 0, 0), picker=True)
        

    #Creates a scatter plot. X is nucleotide position (for selected rrnA operon), Y is read count on specific position.
    if scatter_flag:
        if MA_flag:
            ax1 = ax12
            ax2 = ax22
            
        #Generate 16S scatter plot.
        ax1.set_ylabel('Read Counts', fontsize=14)
        ax1.set_xlabel('Nucleotide Position', fontsize=14)
        ax1.set_title('16S Scatterplot', fontsize=14)

        for data_key, data_nt in data_dic.items():#[data1, data2, data3, data4]:
            ax1.scatter(data_nt.nucl_data_16S, data_nt.y_pos_16S, alpha=0.5, c=data_nt.colour_16S, linewidths=( 0, 0, 0), picker=True, marker = data_nt.symbol)
            ax1.scatter(data_nt.nucl_data_16S, [-1 * data for data in data_nt.y_neg_16S],
                        alpha=0.5, c=data_nt.colour_16S, linewidths=( 0, 0, 0), picker=True, marker = data_nt.symbol)

        #Generate 23S scatter plot.
        ax2.set_ylabel('Read Counts', fontsize=14)
        ax2.set_xlabel('Nucleotide Position', fontsize=14)
        ax2.set_title('23S Scatterplot', fontsize=14)
        
        for data_key, data_nt in data_dic.items():#[data1, data2, data3, data4]:
            ax2.scatter(data_nt.nucl_data_23S, data_nt.y_pos_23S, alpha=0.5, c=data_nt.colour_23S, linewidths=( 0, 0, 0), picker=True, marker = data_nt.symbol)
            ax2.scatter(data_nt.nucl_data_23S, [-1 * data for data in data_nt.y_neg_23S],
                        alpha=0.5, c=data_nt.colour_23S, linewidths=( 0, 0, 0), picker=True, marker = data_nt.symbol)




    fig1.tight_layout()
    fig2.tight_layout()
    fig1.savefig('test5.pdf')
    fig2.savefig('test6.pdf')
    fig1.canvas.mpl_connect('pick_event', onpick_16S)
    fig2.canvas.mpl_connect('pick_event', onpick_23S)

    plt.show()