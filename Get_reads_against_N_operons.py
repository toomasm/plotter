import pandas as pd
import os

def pandas_processer(index, root_ext, N_of_operons):

    genome_index = root_ext[0]
    #Get the file extention input file
    name = genome_index + "_{0}operons.csv".format(N_of_operons)
    #Creates dataframe from csv
    df = pd.read_csv(index)
    
    #Fills empty dataspaces with 0-s
    df = df.fillna(0)

    #Create a new column in dataframe, which contains the number of rRNA operons which read maps to. 
    df['gene_count'] = df['rrnA'] + df['rrnB'] + df['rrnC'] + df['rrnD'] + df['rrnE'] + df['rrnG'] + df['rrnH']
    #Creates a new fataframe which contains only reads which map to only a certain number of rRNA operons
    #selected_df = df.query('gene_count == 7')
    for k in ['Readname', 'rrsA', 'rrsB', 'rrsC', 'rrsD', 'rrsE', 'rrsG', 'rrsH']:
        del df[k]
        
    selected_df = df[df['gene_count'] == N_of_operons]
    #Save the new dataframe
    selected_df.to_csv(name,index=True,header=True)

    return name










