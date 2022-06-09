import seaborn as sns
import matplotlib
import numpy as np
import os
from matplotlib import pyplot as plt
import pandas as pd
import argparse

sns.set_theme() 
#chromosome lengths
chr_lens = {
    'chr1': 248387328,
    'chr2': 242696752,
    'chr3': 201105948,
    'chr4': 193574945,
    'chr5': 182045439,
    'chr6': 172126628,
    'chr7': 160567428,
    'chr8': 146259331,
    'chr9': 150617247,
    'chr10': 134758134,
    'chr11': 135127769,
    'chr12': 133324548,
    'chr13': 113566686,
    'chr14': 101161492,
    'chr15': 99753195,
    'chr16': 96330374,
    'chr17': 84276897,
    'chr18': 80542538,
    'chr19': 61707364,
    'chr20': 66210255,
    'chr21': 45090682,
    'chr22': 51324926,
    'chrX': 154259566,
}
#chromosome names
chr_names = {
    'chr1': 'Chromosome 1',
    'chr2': 'Chromosome 2',
    'chr3': 'Chromosome 3',
    'chr4': 'Chromosome 4',
    'chr5': 'Chromosome 5',
    'chr6': 'Chromosome 6',
    'chr7': 'Chromosome 7',
    'chr8': 'Chromosome 8',
    'chr9': 'Chromosome 9',
    'chr10': 'Chromosome 10',
    'chr11': 'Chromosome 11',
    'chr12': 'Chromosome 12',
    'chr13': 'Chromosome 13',
    'chr14': 'Chromosome 14',
    'chr15': 'Chromosome 15',
    'chr16': 'Chromosome 16',
    'chr17': 'Chromosome 17',
    'chr18': 'Chromosome 18',
    'chr19': 'Chromosome 19',
    'chr20': 'Chromosome 20',
    'chr21': 'Chromosome 21',
    'chr22': 'Chromosome 22',
    'chrX': 'Chromosome X',
}


def main(args):
    #path to the folder with all the .csv files
    input_path = args.data
    #directtory/folder to save all the .png images
    output = args.out
    #create a folder to store all the images/graphs/plots
    if not os.path.isdir(output):
        os.makedirs(output)
    #create bins
    amount_of_bins = 30
    #iterate through all the chromosomes
    for chr_n in chr_lens.keys():
        print(f'==== Processing Plot: {chr_n} ====')
        #path to the csv file
        f = os.path.join(input_path, f'{chr_n}.csv')
        df = pd.read_csv(f)

        #use a scalar bin size otherwise not all the rows are added
        bin_size = chr_lens[chr_n]//amount_of_bins
        bin_array = np.arange(0, chr_lens[chr_n], bin_size)
        
        #Put the edges into bins based on the source start of the edges
        #Method 1: Gives the indices of the intervals in bin arrary
        #df['bin'] = np.digitize(df['Source_Start'], bins=bin_array)
        
        #Method 2: Gives the intervals directly
        df['bin'] = pd.cut(df['Source_Start'] , bins=bin_array, include_lowest=False)
        
        #groups by the bin interval
        df1_grouped = df.groupby('bin')
        
        #create a dictionary for storing the precision & recall for each bin
        bin_dict = {}
        for group_name, df_group in df1_grouped:
            bin_dict[f'{group_name}'] = []
        s = 0

        #iterate through the grouped dataframe
        for group_name, df_group in df1_grouped:
            #calculate the no of TP,FN,TN,FP
            TP = df_group.query('Edge_Type == "TP"').Edge_Type.count()
            FN = df_group.query('Edge_Type == "FN"').Edge_Type.count()
            FP = df_group.query('Edge_Type == "FP"').Edge_Type.count()

            precision = TP / (TP + FN)
 
            #calculate the recall
            recall = TP / (TP + FP)
            
            #add the precision and recall to the dictionary of the interval
            bin_dict[f'{group_name}'] = precision,recall
            
        #start plotting
        #keys = intervals,vals_1 = precision, vals_2 = recall
        keys = []
        vals_1 = []
        vals_2 = []
        keys = list(bin_dict.keys())
        #print([bin_dict[k] for k in keys])
        vals_1 = [float(bin_dict[k][0]) for k in keys]
        vals_2 = [float(bin_dict[k][1]) for k in keys]
            
        #transfer everything to a dataset
        data_preproc = pd.DataFrame({
                'keys': keys, 
                'precision': vals_1,
                'recall': vals_2})
            
        #snslineplot with both precision and recall line
        line_plot = sns.lineplot(x='keys', y='value',hue='variable', data=pd.melt(data_preproc, ['keys']),marker = 'o')
            
        #alter plot aesthetics slightly
        line_plot.set_title(chr_names[chr_n])
        line_plot.set_ylabel("Precision/Recall", fontsize=12)
        line_plot.set_xlabel("Reference Position", fontsize=12)
        plt.setp(line_plot.get_xticklabels(), rotation=90, horizontalalignment='right',fontsize=8)
        plt.ylim(-0.1,1.05)
        #get a .png file to save the plot
        #fig = line_plot.get_figure()
        plt.savefig(f'{output}/{chr_n}.png',bbox_inches="tight")
        plt.clf()

       
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Takes data path')
    parser.add_argument('--data', type=str, default ='excel_path',help='path to the data (folder with chromosome folders)')
    parser.add_argument('--out', type=str,default = 'graph_path', help='path to the output folder')
    args = parser.parse_args()
    main(args)
