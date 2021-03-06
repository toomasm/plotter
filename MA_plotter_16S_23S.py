import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
import pandas as pd
from collections import namedtuple, OrderedDict, Counter, defaultdict
import math
import numpy as np
from Bio import SeqIO
from collections import defaultdict
import hdf5_gen
from progressbar import Bar, ETA, \
    ProgressBar, ReverseBar
import os

LENGTH_OF_MATURE_17S = 1690
#LENGTH_OF_MATURE_23S = 2905
LENGTH_OF_MATURE_23S = 2944



def mature_ribosome_RNA_divider(counter, rRNA):
    """ Takes dictionary containing all rRNA operon positions as keys and read counts in position as values as input.
    Generates a new dictionary representeing a blank 16S or 23S rRNA. Adds all those reads in appropriate position on
    the blank rRNA."""

    #ribosome_dict = {}
    ribosome_dict_pos = defaultdict(float)
    ribosome_dict_neg = defaultdict(float)
    if rRNA == "16S":
        subunit = 0
    elif rRNA == "23S":
        subunit = 1

    for operon in all_operons:

        #Change in positions to compensate the deletions and insertions in rrlA.
        compensate_dif = 0

        if operon in ["rrnD", "rrnG"]:
            operon_type = 'neg'
        else:
            operon_type = 'pos'

        for position, value in counter.items():
            if math.fabs(position) not in range (gene_dic[operon][subunit].start_pos, gene_dic[operon][subunit].end_pos + 1):
                continue

            # rrlA (23S of rrnA operon) has 2 insertion and 1 deletion compared to other 6 rrnA-s.
            #This code is compensating for that while placing reads on our general ribosomal RNA.
            if operon == "rrnA":
                if math.fabs(position) in range (4038692, 4039385):
                    compensate_dif = 1
                elif math.fabs(position) in range (4039391, 4040443+1):
                    compensate_dif = -1

            if operon_type == 'pos':
                new_position = math.fabs(math.fabs(position) - gene_dic[operon][subunit].start_pos) + 1 + compensate_dif
            else:
                new_position = math.fabs(math.fabs(position) - gene_dic[operon][subunit].end_pos) + 1 + compensate_dif

            if (operon_type == 'pos' and position > 0) or (operon_type == 'neg' and position < 0):
                ribosome_dict_pos[new_position] += math.fabs(value)
            else:
                ribosome_dict_neg[new_position] += math.fabs(value)

    return ribosome_dict_pos, ribosome_dict_neg


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



def data_generator(index, standard_colour, ACA_colour, MqsR_colour, symbol):

    global all_operons, gene_dic

    GeneInfo = namedtuple('GeneInfo', ['start_pos', 'end_pos', 'column_label'])
    plot_data = namedtuple('PlottingData', ['nucl_data_16S', 'nucl_data_23S' ,'y_pos_16S', 'y_neg_16S', 'y_pos_23S',
                                            'y_neg_23S', 'colour_16S', 'colour_23S', 'symbol'])

    _16S = namedtuple('_16S', ['start_pos', 'end_pos'])
    _23S = namedtuple('_23S', ['start_pos', 'end_pos'])

    #original 16S and 23S starting and ending position.
    '''gene_dic = OrderedDict([('rrnA', (_16S(4035531, 4037072), _23S(4037519, 4040423))),
                            ('rrnB', (_16S(4166659, 4168200), _23S(4168641, 4171544))),
                            ('rrnC', (_16S(3941808, 3943349), _23S(3943704, 3946607))),
                            ('rrnD', (_16S(3427221, 3428762), _23S(3423880, 3426783))),
                            ('rrnE', (_16S(4208147, 4209688), _23S(4210043, 4212946))),
                            ('rrnG', (_16S(2729616, 2731157), _23S(2726281, 2729184))),
                            ('rrnH', (_16S(223771, 225312), _23S(225759, 228662)))])

    #23S and 16S 20 nucleotides before starting position and 20 nucleotides after ending positions.

    gene_dic = OrderedDict([('rrnB', (_16S(4166639, 4168220), _23S(4168621, 4171564))),
                            ('rrnC', (_16S(3941788, 3943369), _23S(3943684, 3946627))),
                            ('rrnD', (_16S(3427201, 3428782), _23S(3423880, 3426823))),
                            ('rrnE', (_16S(4208127, 4209708), _23S(4210023, 4212966))),
                            ('rrnG', (_16S(2729506, 2731177), _23S(2726281, 2729224))),
                            ('rrnH', (_16S(223751, 225332), _23S(225759, 228682))),
                            ('rrnA', (_16S(4035511, 4037092), _23S(4037499, 4040443)))])'''

    #for 16S the 115 nucleotides before and 33 nucleotides after mature 16S rRNA(first cleavage product by RNaseIII 17S)
    #for 23S 20nucleotides before starting position and 20 nucleotides after ending positions.
    gene_dic = OrderedDict([('rrnB', (_16S(4166544, 4168233), _23S(4168621, 4171564))),
                            ('rrnC', (_16S(3941693, 3943382), _23S(3943684, 3946627))),
                            ('rrnD', (_16S(3427188, 3428877), _23S(3423880, 3426823))),
                            ('rrnE', (_16S(4208032, 4209721), _23S(4210023, 4212966))),
                            ('rrnG', (_16S(2729583, 2731272), _23S(2726281, 2729224))),
                            ('rrnH', (_16S(223656, 225345), _23S(225759, 228682))),
                            ('rrnA', (_16S(4035416, 4037105), _23S(4037499, 4040443)))])

    all_operons = ['rrnB', 'rrnC', 'rrnD', 'rrnE', 'rrnG', 'rrnH', 'rrnA']

    ref_genome_fasta_17S = 'data/17S_new.fasta'
    #ref_genome_fasta_23S = 'rrlA.fasta'
    ref_genome_fasta_23S = 'data/23S_new.fasta'

    dataframe = pd.read_csv(index)

    values_list_pos = []
    values_list_neg = []



    for operon in all_operons:

        number_of_reads = len(dataframe)
        widgets = [Bar('>'), ' ', ETA(), ' ', ReverseBar('<')]
        pbar = ProgressBar(widgets=widgets, maxval=number_of_reads).start()
        i = 0

        for value in dataframe[operon].dropna():

            if '+AC0-' in str(value):
                value = int(float(str(value)[4:]))
                values_list_pos.append(value)
            else:
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
    counter_pos_16S, counter_neg_16S = mature_ribosome_RNA_divider(counter_pos, "16S")
    counter_pos_23S, counter_neg_23S = mature_ribosome_RNA_divider(counter_pos, "23S")

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

    for fasta in SeqIO.parse(ref_genome_fasta_17S, "fasta"):
        pass

    for index in range(1, LENGTH_OF_MATURE_17S + 1):

        if str(fasta.seq[index-1: index+2]) == 'ACA':
            colour_16S.append(ACA_colour)
        elif str(fasta.seq[index-2: index+1]) == 'GCT' or str(fasta.seq[index-2: index+1]) == 'GCA':
            colour_16S.append(MqsR_colour)
        else:
            colour_16S.append(standard_colour)
        nucl_data_16S.append(index-115)
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



def make_plot_all_reads(_input_, scatter_flag, delta_scatter_flag, MA_flag, three_prime_flag):
    """ Creates a MA plot or scatter plot against chosen rrnA operon from two CSV index files (index, positions, mapped against how many).
        Because my system has memory limitations with pandas I strip down the indexes before plotting them here. All reads data
        I want to show is present in the CSV file."""
    dataset_names = []
    processing_names = []
    primes = []
    for ident in _input_:
        ident = os.path.basename(ident)
        ident = ident.split('_')
        name = ident[0]
        process = ident[1]
        prime = ident[2]
        dataset_names.append(name)
        processing_names.append(process)
        primes.append(prime)
    print primes
    print dataset_names
    print processing_names


    try:
        #looks how many input files have been given
        input_number = len(_input_)
        data_dic = hdf5_gen.get_hdf_data(primes, dataset_names, processing_names, _input_)
       

    except (KeyError, IOError):
        if three_prime_flag:
            symbol_1 = 'o'
            symbol_2 = 'D'

        else:
            symbol_1 = 'o'
            symbol_2 = 'o'

        hdf_dic = defaultdict(lambda : defaultdict(dict))

        #Creates data necessary for plotting for every inputfile.
        hdf_dic[primes[0]][dataset_names[0]][processing_names[0]] = data_generator(_input_[0], 'b', 'r', 'm', symbol_1)
        if input_number >= 2:
            hdf_dic[primes[1]][dataset_names[1]][processing_names[1]] = data_generator(_input_[1], 'b', 'r', 'm', symbol_2)

        if input_number >= 3:
            hdf_dic[primes[2]][dataset_names[2]][processing_names[2]] = data_generator(_input_[2], 'b', 'r', 'm', symbol_1)

        if input_number == 4:
            hdf_dic[primes[3]][dataset_names[3]][processing_names[3]] = data_generator(_input_[3], 'b', 'r', 'm', symbol_2)


        hdf5_gen.add_to_hdf(hdf_dic, 'hdf5/plot_data.hdf')
        data_dic = hdf5_gen.get_hdf_data(primes, dataset_names,processing_names,_input_)

    def onpick_16S(event):

        #Gets the index of datapoint (index of X and Y)
        ind = event.ind

        #Retrieves information from lists based on index of the datapoint.
        #The information to be displayed for datapoint has the same index as datapint.

        nucleotide_pos = np.take(data_dic['data1']['nucl_data_16S'], ind)
        sample_1_pos_read_count = np.take(data_dic['data1']['y_pos_16S'], ind)
        sample_2_pos_read_count = np.take(data_dic['data2']['y_pos_16S'], ind)
        sample_1_neg_read_count = np.take(data_dic['data1']['y_neg_16S'], ind)
        sample_2_neg_read_count = np.take(data_dic['data2']['y_neg_16S'], ind)

        if input_number >= 3:
            sample_3_pos_read_count = np.take(data_dic['data3']['y_pos_16S'], ind)
            sample_3_neg_read_count = np.take(data_dic['data3']['y_neg_16S'], ind)
        if input_number == 4:
            sample_4_pos_read_count = np.take(data_dic['data4']['y_pos_16S'], ind)
            sample_4_neg_read_count = np.take(data_dic['data4']['y_neg_16S'], ind)

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

        nucleotide_pos = np.take(data_dic['data1']['nucl_data_23S'], ind)
        sample_1_pos_read_count = np.take(data_dic['data1']['y_pos_23S'], ind)
        sample_2_pos_read_count = np.take(data_dic['data2']['y_pos_23S'], ind)
        sample_1_neg_read_count = np.take(data_dic['data1']['y_neg_23S'], ind)
        sample_2_neg_read_count = np.take(data_dic['data2']['y_neg_23S'], ind)

        if input_number >= 3:
            sample_3_pos_read_count = np.take(data_dic['data3']['y_pos_23S'], ind)
            sample_3_neg_read_count = np.take(data_dic['data3']['y_neg_23S'], ind)
        if input_number == 4:
            sample_4_pos_read_count = np.take(data_dic['data4']['y_pos_23S'], ind)
            sample_4_neg_read_count = np.take(data_dic['data4']['y_neg_23S'], ind)

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
    if MA_flag and (scatter_flag or delta_scatter_flag):
        fig1, (ax11, ax12) = plt.subplots(2)
        fig2, (ax21, ax22) = plt.subplots(2)
    else:
        fig1, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots()

    if MA_flag:
        if scatter_flag or delta_scatter_flag:
            ax1 = ax11
            ax2 = ax21

        #Generates data for MA plot.
        MA_X_16S = [(math.log(float(y1), 2) + math.log(float(y2), 2))/2 for y1, y2 in zip(data_dic['data1']['y_pos_16S'], data_dic['data2']['y_pos_16S'])]
        MA_Y_16S = [math.log((float(y1)/y2), 2) for y1, y2 in zip(data_dic['data1']['y_pos_16S'], data_dic['data2']['y_pos_16S'])]

        MA_X_23S = [(math.log(float(y1), 2) + math.log(float(y2), 2))/2 for y1, y2 in zip(data_dic['data1']['y_pos_23S'], data_dic['data2']['y_pos_23S'])]
        MA_Y_23S = [math.log((float(y1)/y2), 2) for y1, y2 in zip(data_dic['data1']['y_pos_23S'], data_dic['data2']['y_pos_23S'])]


        ax1.set_ylabel('M', fontsize=20)
        ax1.set_xlabel('A', fontsize=20)
        ax1.set_title('16S MA plot', fontsize=14)
        ax1.scatter(MA_X_16S, MA_Y_16S, alpha=0.5, c=data_dic["data1"]['colour_16S'], linewidths=( 0, 0, 0), picker=True)

        ax2.set_ylabel('M', fontsize=20)
        ax2.set_xlabel('A', fontsize=20)
        ax2.set_title('23S MA plot', fontsize=14)
        ax2.scatter(MA_X_23S, MA_Y_23S, alpha=0.5, c=data_dic["data1"]['colour_23S'], linewidths=( 0, 0, 0), picker=True)


    #Creates a scatter plot. X is nucleotide position, Y is 5 prime read end count or 5 prime read counts relative percentage
    #compared against 5 prime end of mature rRNA-s (16S or 23S).
    if scatter_flag:
        if MA_flag:
            ax1 = ax12
            ax2 = ax22

        #Set title and x-label for 16S and 23S scatterplot
        ax1.set_xlabel('Nucleotide Position', fontsize=14)
        ax1.set_title('16S Scatterplot', fontsize=14)

        ax2.set_xlabel('Nucleotide Position', fontsize=14)
        ax2.set_title('23S Scatterplot', fontsize=14)

        if scatter_flag == "percent":

            ax1.set_ylabel('Relative percentage of reads', fontsize=14)
            ax2.set_ylabel('Relative percentage of reads', fontsize=14)

            hundred_percent_16S = 0
            hundred_percent_23S = 0

            for data_key, data_nt in data_dic.items():#[data1, data2, data3, data4]:

                for i in range(-19,22):
                
                    index = data_dic['data1']['nucl_data_16S'].index(float(i))
                    hundred_percent_16S += data_nt['y_pos_16S'][index]
                    hundred_percent_23S += data_nt['y_pos_23S'][index]
                    
                heights_16S = [(x / hundred_percent_16S * 100) for x in data_nt['y_pos_16S']]
                heights_23S = [(x / hundred_percent_23S * 100) for x in data_nt['y_pos_23S']]
                
                ax1.scatter(data_nt['nucl_data_16S'], heights_16S, alpha=0.5, c=data_nt['colour_16S'], linewidths=( 0, 0, 0),
                            picker=True, marker = data_nt['symbol'])
                ax2.scatter(data_nt['nucl_data_23S'], heights_23S, alpha=0.5, c=data_nt['colour_23S'], linewidths=( 0, 0, 0),
                            picker=True, marker = data_nt['symbol'])
                
                """df_16S_pos_none = pd.read_csv("16S_count_pos_none.csv", index_col=0, header= [0,1])
                df_23S_pos_none = pd.read_csv("23S_count_pos_none.csv", index_col=0, header= [0,1]) 
                df_16S_pos_PNK = pd.read_csv("16S_count_pos_PNK.csv", index_col=0, header= [0,1])
                df_23S_pos_PNK = pd.read_csv("23S_count_pos_PNK.csv", index_col=0, header= [0,1])

                #a = data_dic['data1']['y_pos_16S']
                #b = data_dic['data1']['y_pos_23S']
                df_16S_pos_none["MG_log", "count"] = data_dic['data2']['y_pos_16S']
                df_23S_pos_none["MG_log", "count"] = data_dic['data2']['y_pos_23S']
                df_16S_pos_PNK["MG_log", "count"] = data_dic['data1']['y_pos_16S']
                df_23S_pos_PNK["MG_log", "count"] = data_dic['data1']['y_pos_23S']

                df_16S_pos_none.to_csv("16S_count_pos_none.csv", index=True ,header=True)
                df_23S_pos_none.to_csv("23S_count_pos_none.csv", index=True ,header=True)
                df_16S_pos_PNK.to_csv("16S_count_pos_PNK.csv", index=True ,header=True)
                df_23S_pos_PNK.to_csv("23S_count_pos_PNK.csv", index=True ,header=True)"""
              
        elif scatter_flag == "count":
            
            #Set y axis label for scatterplot with read counts.
            ax1.set_ylabel('Read Counts', fontsize=14)
            ax2.set_ylabel('Read Counts', fontsize=14)

            #Makes scatterplots for each datasets 16S and 23S rRNA.
            for data_key, data_nt in data_dic.items():#[data1, data2, data3, data4]:

                ax1.scatter(data_nt['nucl_data_16S'], data_nt['y_pos_16S'], alpha=0.5, c=data_nt['colour_23S'], linewidths=( 0, 0, 0),
                            picker=True, marker = data_nt['symbol'])
                ax1.scatter(data_nt['nucl_data_16S'], [-1 * data for data in data_nt['y_neg_16S']],
                            alpha=0.5, c=data_nt['colour_23S'], linewidths=( 0, 0, 0), picker=True, marker = data_nt['symbol'])

                ax2.scatter(data_nt['nucl_data_23S'], data_nt['y_pos_23S'], alpha=0.5, c=data_nt['colour_23S'], linewidths=( 0, 0, 0),
                            picker=True, marker = data_nt['symbol'])
                ax2.scatter(data_nt['nucl_data_23S'], [-1 * data for data in data_nt['y_neg_23S']],
                            alpha=0.5, c=data_nt['colour_23S'], linewidths=( 0, 0, 0), picker=True, marker = data_nt['symbol'])
        
            '''df_16S_pos_none = pd.read_csv("tables/17S_count_pos_TAP.csv", index_col=0, header= [0,1])
            df_23S_pos_none = pd.read_csv("tables/23S_count_pos_TAP_2.csv", index_col=0, header= [0,1]) 
            df_16S_pos_PNK = pd.read_csv("tables/17S_count_neg_TAP.csv", index_col=0, header= [0,1])
            df_23S_pos_PNK = pd.read_csv("tables/23S_count_neg_TAP_2.csv", index_col=0, header= [0,1])

            #a = data_dic['data1']['y_pos_16S']
            #b = data_dic['data1']['y_pos_23S']
            df_16S_pos_none["MG_log", "count"] = data_dic['data3']['y_pos_16S']
            df_23S_pos_none["MG_log", "count"] = data_dic['data3']['y_pos_23S']
            df_16S_pos_PNK["MG_log", "count"] = data_dic['data3']['y_neg_16S']
            df_23S_pos_PNK["MG_log", "count"] = data_dic['data3']['y_neg_23S']

            df_16S_pos_none.to_csv("tables/17S_count_pos_TAP.csv", index=True ,header=True)
            df_23S_pos_none.to_csv("tables/23S_count_pos_TAP_2.csv", index=True ,header=True)
            df_16S_pos_PNK.to_csv("tables/17S_count_neg_TAP.csv", index=True ,header=True)
            df_23S_pos_PNK.to_csv("tables/23S_count_neg_TAP_2.csv", index=True ,header=True)'''

            count_reads_1_16S = 0
            count_reads_2_16S = 0
            count_reads_1_23S = 0
            count_reads_2_23S = 0
        
        for count_1_16S in data_dic['data1']['y_pos_16S']:#[0:115 + 1 + 20 + 1]:#, data_dic['data2']['y_pos_16S']:
            count_reads_1_16S += count_1_16S
            #count_reads_2 += count_2
            
        for count_2_16S in data_dic['data2']['y_pos_16S']:#[0:115 + 1 + 20 + 1]:#, data_dic['data2']['y_pos_16S']:
            count_reads_2_16S += count_2_16S
            
        for count_1_23S in data_dic['data1']['y_pos_23S']:#[0:42]:#, data_dic['data2']['y_pos_16S']:
            count_reads_1_23S += count_1_23S
            #count_reads_2 += count_2
            
        for count_2_23S in data_dic['data2']['y_pos_23S']:#[0:42]:#, data_dic['data2']['y_pos_16S']:
            count_reads_2_23S += count_2_23S
            
        print count_reads_1_16S
        print count_reads_2_16S
        print count_reads_1_23S
        print count_reads_2_23S
        
    if delta_scatter_flag:
        if MA_flag:
           ax1 = ax12
           ax2 = ax22

        data_nt_pos_16S = [x1 - x2 for x1, x2 in zip(data_dic['data1']['y_pos_16S'], data_dic['data2']['y_pos_16S'])]
        data_nt_pos_23S = [x1 - x2 for x1, x2 in zip(data_dic['data1']['y_pos_23S'], data_dic['data2']['y_pos_23S'])]
        data_nt_neg_16S = [x1 - x2 for x1, x2 in zip(data_dic['data1']['y_neg_16S'], data_dic['data2']['y_neg_16S'])]
        data_nt_neg_23S = [x1 - x2 for x1, x2 in zip(data_dic['data1']['y_neg_23S'], data_dic['data2']['y_neg_23S'])]

        #dict1 = []
        #[dict1.append((x, math.fabs(y))) for x, y in zip(data_dic['data1']['nucl_data_16S'], data_nt_pos_16S)]
        #del dict1[0]
        #print sorted(dict1, key=lambda pos: pos[1],  reverse = True)[:15]

        """df_16S_pos = pd.read_csv("16S_table_pos.csv", index_col=0)
        df_23S_pos = pd.read_csv("23S_table_pos.csv", index_col=0) 
        df_16S_neg = pd.read_csv("16S_table_neg.csv", index_col=0) 
        df_23S_neg = pd.read_csv("23S_table_neg.csv", index_col=0)
       
        df_16S_pos["d3_PNK/d3_none"] = data_nt_pos_16S
        df_23S_pos["d3_PNK/d3_none"] = data_nt_pos_23S
        df_16S_neg["d3_PNK/d3_none"] = data_nt_neg_16S
        df_23S_neg["d3_PNK/d3_none"] = data_nt_neg_23S

        df_16S_pos.to_csv("16S_table_pos.csv", index=True ,header=True)
        df_23S_pos.to_csv("23S_table_pos.csv", index=True ,header=True)
        df_16S_neg.to_csv("16S_table_neg.csv", index=True ,header=True)
        df_23S_neg.to_csv("23S_table_neg.csv", index=True ,header=True)"""
        
        """df_16S_pos = pd.read_csv("16S_pos_none.csv", index_col=0, header= [0,1])
        df_23S_pos = pd.read_csv("23S_pos_none.csv", index_col=0, header= [0,1]) 
        df_16S_neg = pd.read_csv("16S_neg_none.csv", index_col=0, header= [0,1]) 
        df_23S_neg = pd.read_csv("23S_neg_none.csv", index_col=0, header= [0,1])
       
        df_16S_pos["d3_none/MG_log_none", "delta_count"] = data_nt_pos_16S
        df_23S_pos["d3_none/MG_log_none", "delta_count"] = data_nt_pos_23S
        df_16S_neg["d3_none/MG_log_none", "delta_count"] = data_nt_neg_16S
        df_23S_neg["d3_none/MG_log_none", "delta_count"] = data_nt_neg_23S

        df_16S_pos.to_csv("16S_pos_none.csv", index=True ,header=True)
        df_23S_pos.to_csv("23S_pos_none.csv", index=True ,header=True)
        df_16S_neg.to_csv("16S_neg_none.csv", index=True ,header=True)
        df_23S_neg.to_csv("23S_neg_none.csv", index=True ,header=True)"""


        #Generate 16S scatter plot.
        ax1.set_ylabel('Read Count Difference', fontsize=14)
        ax1.set_xlabel('Nucleotide Position', fontsize=14)
        ax1.set_title('16S Scatterplot', fontsize=14)


        ax1.scatter(data_dic['data1']['nucl_data_16S'], data_nt_pos_16S, alpha=0.5, c=data_dic['data1']['colour_16S'], linewidths=( 0, 0, 0), picker=True,
                    marker = data_dic['data1']['symbol'])
        ax1.scatter(data_dic['data1']['nucl_data_16S'], [-1 * data for data in data_nt_neg_16S],
                    alpha=0.5, c=data_dic['data1']['colour_16S'], linewidths=( 0, 0, 0), picker=True, marker = data_dic['data1']['symbol'])

        #Generate 23S scatter plot.
        ax2.set_ylabel('Read Count Difference', fontsize=14)
        ax2.set_xlabel('Nucleotide Position', fontsize=14)
        ax2.set_title('23S Scatterplot', fontsize=14)

        ax2.scatter(data_dic['data1']['nucl_data_23S'], data_nt_pos_23S, alpha=0.5, c=data_dic['data1']['colour_23S'], linewidths=( 0, 0, 0), picker=True,
                    marker = data_dic['data1']['symbol'])
        ax2.scatter(data_dic['data1']['nucl_data_23S'], [-1 * data for data in data_nt_neg_23S],
                    alpha=0.5, c=data_dic['data1']['colour_23S'], linewidths=( 0, 0, 0), picker=True, marker = data_dic['data1']['symbol'])


    fig1.tight_layout()
    fig2.tight_layout()
    fig1.savefig('test5.pdf')
    fig2.savefig('test6.pdf')
    fig1.canvas.mpl_connect('pick_event', onpick_16S)
    fig2.canvas.mpl_connect('pick_event', onpick_23S)

    plt.show()

#comment comment
