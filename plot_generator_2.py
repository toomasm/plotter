import textwrap
import os
import argparse
from collections import namedtuple
from MA_plotter_2 import by_strand_sorter, data_generator, make_plot
from MA_plotter_16S_23S import make_plot_all_reads
from Get_reads_against_N_operons import pandas_processer
#from plotter import read_trimmer, three_prime_checker
#from align_and_index import perform_alignment, sam_to_bam
#from index_generator import generate_index

def make_argument_parser():
    '''Returns argument parser for this script.
    '''
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                             description=textwrap.dedent('''
                             This script takes indexed CSV files (index, readname, rRNA operons, positions) as an input and creates scatter plots
                             (X = position, Y = read count) MA scatter plots. You can choose against which MG1655 E. coli operon (Genbank version 913.3)
                             you want the reads to be plotted.
                                                         
                             On scatter plot you can plot both negative and positive strand reads (negative strand will be on -Y axis).
                             You can specify if you want to show reads that map to all operons or only to some.
                             '''))

    operon_args = parser.add_argument_group('Specify reads chosen for plotting based on their occurende in rRNA operons')

    operon_args.add_argument('-n','--mapped-to-N-operons',
                              default=None, type=int, dest='mapped_to_N_operons',
                              help='Sets how many rRNA operons the reads should map to. Input is CSV (index, readname, rRNA operons, positions)\
                                    and output is a new CSV file that only contains reads which fulfill  the criteria')

    operon_args.add_argument('-o','--show-on-operon',
                              default='rrnA', dest = 'show_on_operon',
                              help='Specify the rrnA operon against which the data will be plotted. Default is rrnA.\
                                    If you specifify that reads have to map against less than 7 operons, all the reads will have to map\
                                    against this one.\n\
                                                     \n\
                                    Choises are: rrnA\n\
                                                 rrnB\n\
                                                 rrnC\n\
                                                 rrnD\n\
                                                 rrnE\n\
                                                 rrnG\n\
                                                 rrnH')

    plotting_args = parser.add_argument_group('Plotting options.')
    
    plotting_args.add_argument('-s', '--scatter-plot', action = 'store_true',
                                dest='scatter_plot', default=None, 
                                help='Makes a scatter plot from CSV file.')
    
    plotting_args.add_argument('-MA', '--MA-plot', action = 'store_true',
                                dest='MA_plot', default=None, 
                                help='Makes a MA scatter plot from 2 CSV files. Two files required for input. Use -i2')

    plotting_args.add_argument('-d', '--delta-scatter-plot', action = 'store_true',
                                dest='delta_scatter_plot', default=None, 
                                help='Makes a scatter plot from CSV file.')
    
    plotting_args.add_argument('-3', '--three_prime', action = 'store_true',
                                dest='three_prime', default=False, 
                                help='Flag that sets parameters for scatter plot to show three prime seq data. Poly \
                                     A tail problem has led us to generate ok and shady reads for each 3 prime treatment.\n\n\
                                     Input has to be always 4 files if this flag is set. The program assumas that you enter them in \
                                     order of sample 1 ok, sample 1 shady, sample 2 ok, sample 2 shady.')

    plotting_args.add_argument('-a', '--plot_all_reads', action = 'store_true',
                                dest='plot_all_reads', default=False, 
                                help='Plots all reads that have been mapped to any rrnA operon against 16S and 23S rRNA')

    input_args = parser.add_argument_group('Data input options.')
    
    input_args.add_argument('-i', '--input-file', 
                                dest='input_filename', default=None, required = True, nargs = '+',
                                help='Input filename. This argument takes in upto 4 CSV files. Used for plotting two to four datasets upon one another \
                                      or one dataset by itself.\
                                      Required argument')

    
    return parser


    
def workdir_and_names(input_name):
    """ Gets the basename, root extention and file extention of input file. """

    input_data = namedtuple('input_data', ['base_name' ,'root_ext', 'file_type'])

    # Get the root file name for the input file.
    base_name = os.path.basename(input_name)

    # Get the root and extention of the base_name.
    root_ext = os.path.splitext(base_name)

    #Get the extention of input file.
    input_filetype = root_ext[-1]

    return input_data(base_name, root_ext, input_filetype)


    
def process_arguments(args):

    #Changes working directory to first input files directory.
    #abspath = os.path.abspath(args.input_filename[0])
    #dir_name = os.path.dirname(abspath)
    #os.chdir(dir_name)
    
    #Checks that there are no more inputs than 4.
    if len(args.input_filename) > 4:
        raise UserWarning ('Maximum number of inputs is 4!')

    #Count for generating differenet input variable names.
    count = 0
    #List of variable names that represent files paramameters saved in a namedtuple.
    input_list = []
    #Iterates over input list and generates a filename for each input.
    for inputs in args.input_filename:
        count += 1
        name_of_input = {}
        name = 'input' + str(count)
        #Generates basename, root extention and gets filetype for input file.
        name_of_input[name] = workdir_and_names(inputs)
        #Adds every generated input file data holding tuple to a list.
        input_list.append(name_of_input[name])
        #Checks if input file is  csv file.
        if name_of_input[name].file_type != '.csv':
            raise UserWarning ('Inputs must be CSV files.')
            
    #Generates a directory for plots if not present in working directory or is not the working directory.
    if not os.path.exists('plots'):
        os.makedirs('plots')

    #Generates new CSV files that contain only reads mapping to specific number of operons.
    if args.mapped_to_N_operons:
        #List of newly generated CSV filenames that are optimized for ploting.
        CSV_N = []
        #Iteretaes over the inputs data tuple list, generates a plotting optimised csv file for each inputfile.
        #Output is the path to that file. Adds it to list.
        for inputs in input_list:
            CSV_N.append(pandas_processer(inputs.base_name, inputs.root_ext, args.mapped_to_N_operons))
    else:
        #If inputfiles are already ok for plotting makes the inputfile list the list to contain plotting data. 
        CSV_N = args.input_filename


    #When one plotting argument is set, check if operon has been set to map against.
    if args.scatter_plot or args.MA_plot or args.delta_scatter_plot:
            
        if args.input_filename and args.plot_all_reads:
            make_plot_all_reads(CSV_N, args.scatter_plot, args.delta_scatter_plot, args.MA_plot, args.three_prime)
        elif args.input_filename:
            
            make_plot(CSV_N, args.show_on_operon, args.scatter_plot, args.delta_scatter_plot, args.MA_plot, args.three_prime)

            
if __name__ == '__main__':

    p = make_argument_parser()
    args = p.parse_args()
    process_arguments(args)
    