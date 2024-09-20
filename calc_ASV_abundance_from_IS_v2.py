# Calculate the absolute ASV abundance from the internal standards.
# Previous script only calculated BP, making more versatile to calculate ASV abundance from any internal standard.
#
# Format: calc_ASV_abundance_from_IS_v2.py IS, where IS can be 'BP','DR', or 'TT'
#
# Lexi Jones-Kellett
# Date created: 08/11/23
# Last edited: 09/20/24

import csv,sys,json
import numpy as np

IS = str(sys.argv[1]) # user input is one of the internal standards

if IS == 'BP':
    genome_len = 6244976 #basepairs
    rrn = 5 #[16S rRNA copy#/cell]
    IS_OTUs = ['0a1e7e4b25a59be69931c5d7f92751f5','f40b1be49d3bca5b8fabdd944abb31bf','2029a1010d7bebac2d09361c275f9fda']
elif IS == 'DR':
    genome_len = 3279485
    rrn = 3
    IS_OTUs = ['6a5fcf5f0ca1f18bca2297194442a6d7']
elif IS == 'TT':
    genome_len = 2143708
    rrn = 2
    IS_OTUs = ['9aa3ebacc998945a0cd514ca909e5231','5b0d64b13238ee1991c15a9913bec9bc']
else:
    print('Error')

    
# config file contains local directory paths for "data_dir" and "ASV_count_dir"
with open("conf.json") as json_conf : 
    config = json.load(json_conf)
    
# Get the batch number for each of the samples
sample_batches = {}
with open(config["data_dir"] + 'G4_IS_sample_batches.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader: 
        if i == 0: # header
            pass
        else:
            sample_batches[row[1]] = int(row[0])
        i += 1

# Get the concentrations of the internal standards for each batch (from DNA quantification)
IS_quants = {}
with open(config["data_dir"] + 'IS_quants.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    i = 0
    for row in csv_reader:
        if i == 0:
            pass # skip header
        else:
            if IS == 'BP':
                IS_quants[int(row[0])] = float(row[1]) # BP quant is column index 1
            elif IS == 'DR':
                IS_quants[int(row[0])] = float(row[2]) 
            elif IS == 'TT':
                if int(row[0]) == 2: # exclude batch 2 due to pipetting error, make NAN a multiplier
                    IS_quants[int(row[0])] = np.nan
                else:
                    IS_quants[int(row[0])] = float(row[3])
        i += 1

def calc_IS_Cs(genome_len,rrn,conc):
    """
    Inputs
        genome_len: genome length [bp]
        rrn: copy number of 16S rRNA per cell
        conc: concentration of DNA from quantification, batch specific [ng/muL]
    Output
        Cs: total number of 16S rRNA copies spiked into the samples
    """
    AgC = 6.022*(10**23) #Avogadro's constant [copies/mol]
    bp_weight = 650 #[g/(mol*bp)]
    CF = 10**9 #[ng/g]
    vol_IS = 20 #[muL]
    return (conc*vol_IS*AgC*rrn)/(genome_len*CF*bp_weight)

# Calculate the total number of 16S rRNA copies spiked into each sample
Cs = {}
for key in IS_quants:
    Cs[key] = calc_IS_Cs(genome_len,rrn,IS_quants[key])

# Merged tables
ASV_count_dir = '/Users/lexijones/Dropbox (MIT)/Grad_School/Research/G4_consolidated/data/HighCoverage_221107/merged_16S_18S_tables/'
ASV_count_file = '221118-1309_LexiGradients-HighCov_2.09-fold-18S-correction_normalized_sequence_counts.tsv'
ASV_IS_ABS_file = '221118-1309_LexiGradients-HighCov_2.09-fold-18S-correction_normalized_sequence_counts_abs_ASV_abundance_IS_%s.tsv'%(IS)

IS_counts = []
with open(config["ASV_count_dir"] + ASV_count_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    for row in csv_reader:
        OTU = row[0]
        if OTU in IS_OTUs:
            IS_counts.append([float(r) for r in row[1:-1]]) # exclude OTU and taxonomy
IS_counts_total = np.nansum(IS_counts,axis=0)

###### CREATE NEW TABLE W/ ABS ASV ESTIMATE ####### 

# includes BP,DR,TT to be excluded from count data to save
all_IS_OTUs = ['0a1e7e4b25a59be69931c5d7f92751f5','f40b1be49d3bca5b8fabdd944abb31bf','2029a1010d7bebac2d09361c275f9fda',
           '6a5fcf5f0ca1f18bca2297194442a6d7','9aa3ebacc998945a0cd514ca909e5231','5b0d64b13238ee1991c15a9913bec9bc'] 

f = open(config["ASV_count_dir"] + ASV_IS_ABS_file, 'w') # new table file name
writer = csv.writer(f)

with open(ASV_count_dir + ASV_count_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter='\t')
    i = 0
    for row in csv_reader:
        cropped_row = row[1:-1] # excludes OTU and taxonomy
        if (i == 0):
            writer.writerow(row)
            header = cropped_row
        else:
            OTU = row[0]
            if OTU in all_IS_OTUs: # don't write internal standards into new table
                pass
            else:
                new_row = [OTU]
                for j in np.arange(0,len(cropped_row)): #calculate absolute abundance
                    key = header[j]
                    if key not in sample_batches: #PCR blank samples
                        new_row.append(np.nan)
                    else:
                        if float(IS_counts_total[j]) == 0:
                            new_row.append(np.nan)
                        else:
                            if key == '17C': #950 mL (sample where a bit of volume spilled out)
                                new_row.append((float(cropped_row[j])*float(Cs[sample_batches[key]]))/(float(IS_counts_total[j])*950))
                            else: #1000 mL
                                new_row.append((float(cropped_row[j])*float(Cs[sample_batches[key]]))/(float(IS_counts_total[j])*1000))
                new_row.append(row[-1])
                writer.writerow(new_row)
        i += 1 
  
f.close()

