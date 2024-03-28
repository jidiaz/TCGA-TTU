#!/usr/bin/env python
import os
import pandas as pd
import datetime
from functools import reduce
import argparse
import json

# Dynamically build the path to the config.json file
script_dir = os.path.dirname(os.path.abspath(__file__))  # Gets the directory where the script is located
config_path = os.path.join(script_dir, 'config.json')  # Builds the path to config.json in the same directory
# Load JSON configuration as default settings
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

# Set up argparse with the ability to override JSON settings

parser = argparse.ArgumentParser(description='Concatenates TCGA transcriptome raw counts files')
parser.add_argument("-p", default=config['path_to_TCGA_files'], action='store', help='path to the TCGA files')
parser.add_argument("-d", default=config['clinical_data_path'], type=argparse.FileType('r'), help='clinical data')
parser.add_argument("-r", default=config['sample_sheet_path'], type=argparse.FileType('r'), help='sample sheet from TCGA')
parser.add_argument("-n", action='store_true', default=config['rename_files'], help='A boolean flag to point if the file needs to be renamed, changing the TCGA code for the sample name')
parser.add_argument("-o", default=config['output_path'], action='store', help='path to save files')

args = parser.parse_args()
clinical=pd.read_csv(args.d,sep='\t')
sampleinfo=args.r
path=args.p
folders=args.o
name_change=args.n

id=[]
change=[]

if name_change:
	for filename in os.listdir(path):
		filename = filename.strip('\n')
		id.append(filename)

	for secu in sampleinfo:
		change.append(secu)

	for iname in id:
		siname = iname.split('.rna')[0]
		for fname in change:
			if siname in fname:
				z = fname.split('\t')[6]
				j = z.split(",")
				os.rename(path+iname, path+j[0])

else:
    print("No name change")


#Date used to name output files
current_datetime = datetime.datetime.now()
date_format = "%Y%m%d"
formatted_date = current_datetime.strftime(date_format)
date_in_name = formatted_date[2:]


#extract the name of the project
poject_in_name=clinical['project_id'][0]	

#create the folder where save the files
name_in_output= date_in_name+"_"+poject_in_name
os.makedirs(os.path.join(folders + name_in_output))


#create two tables: TPM+Clinical and raw caount
tpm_clinical={}
raw_count={}
# Initialize an empty list to keep track of skipped files
skipped_files = []

for filename in os.listdir(path):
    # Skip .DS_Store files
    if filename == '.DS_Store':
        skipped_files.append(filename)
        continue

    try:
        # Load each file
        filename = filename.strip('\n')
        tcga_file = pd.read_csv(path + filename, header=None)
        
        # Ensure the data is treated as strings for processing
        tcga_file[0] = tcga_file[0].astype(str)
        
        # Split the first column by tab, expanding into multiple columns
        split_columns = tcga_file[0].str.split('\t', expand=True)
        headers = ['gene_id', 'gene_name', 'gene_type', 'unstranded', 'stranded_first', 'stranded_second', 'tpm_unstranded', 'fpkm_unstranded', 'fpkm_uq_unstranded']
        #print(filename)
        #print("Number of columns after split:", split_columns.shape[1])
        #print("Number of headers:", len(headers))
        if len(headers) != split_columns.shape[1]:
            print("Mismatch in the number of columns and headers.")
            raise ValueError("The number of split columns does not match the number of headers.")

        # Assigning split columns to tcga_file with headers
        tcga_file = split_columns
        tcga_file.columns = headers

        tcga_file = tcga_file.iloc[2:]

        # Processing data specific to your needs
        tcga_file = tcga_file[['gene_id', 'unstranded', 'tpm_unstranded']]
        tcga_file = tcga_file.rename(columns={'unstranded': filename})
        tpm = tcga_file[['gene_id', 'tpm_unstranded']]
        tpm = tpm.rename(columns={'tpm_unstranded': filename})
        tpm_clinical[filename] = tpm
        raw_count[filename] = tcga_file[['gene_id', filename]]

    except Exception as e:
        print(f"Error processing {filename}: {e}")
        skipped_files.append(filename)

# Check if there are any skipped files and write them to a text file
if skipped_files:
    # Define the path for the skipped files text file within the output directory
    skipped_files_path = os.path.join(folders + name_in_output, 'skipped_files.txt')
    
    # Write the names of the skipped files to the text file
    with open(skipped_files_path, 'w') as f:
        for file in skipped_files:
            f.write(file + '\n')
    
    print(f"Skipped files list saved to {skipped_files_path}")


print("initializing merging process, this may take few minutes...")
#Merge files
raw_count_table = reduce(lambda left, right: pd.merge(left, right, on=['gene_id'], how='left'), raw_count.values())
tpm_count_table = reduce(lambda left, right: pd.merge(left, right, on=['gene_id'], how='left'), tpm_clinical.values())
	


#Function to recover genes of interest 
def tpmoverview (dataframe):
	mage_ensmbl = ['ENSG00000198681','ENSG00000268606','ENSG00000183305','ENSG00000221867','ENSG00000147381',
				   'ENSG00000242520','ENSG00000197172','ENSG00000224732'
		,'ENSG00000156009','ENSG00000123584','ENSG00000267978','ENSG00000124260','ENSG00000185247',
				   'ENSG00000213401','ENSG00000214107','ENSG00000099399','ENSG00000198798'
		,'ENSG00000120289','ENSG00000188408','ENSG00000176746','ENSG00000177689','ENSG00000189023',
				   'ENSG00000182798','ENSG00000176774'
		,'ENSG00000155495','ENSG00000046774','ENSG00000165509','ENSG00000179222','ENSG00000102316',
				   'ENSG00000154545','ENSG00000187243','ENSG00000198934'
		,'ENSG00000186675','ENSG00000177383','ENSG00000187601','ENSG00000254585','ENSG00000268916',
		'ENSG00000198930','ENSG00000268902','ENSG00000242599','ENSG00000130726']

	list = dataframe.gene_id.tolist()
	for x in list:
		z = x.split(".")[0]
		dataframe.gene_id = dataframe.gene_id.replace(x, z)

	dataframe = dataframe.loc[dataframe['gene_id'].isin(mage_ensmbl)]


	return(dataframe)

#Select only the genes we are interested in
tpm_merge_ensmblcode = tpmoverview(tpm_count_table)


#
tpm_merge_ensmblcode = tpm_merge_ensmblcode.set_index(['gene_id'])
tpm_merge_ensmblcode = tpm_merge_ensmblcode.T
tpm_merge_ensmblcode = tpm_merge_ensmblcode.reset_index()
codes = tpm_merge_ensmblcode['index']

new_id = []
for x in codes:
    z = x.split('_normal')
    y = z[0].split('-')
    change = y[0] + '-' + y[1] + '-' + y[2]
    new_id.append(change)

#modify the sample ID (delete the number of the replicate) to match with the ID in the clinical file
tpm_merge_ensmblcode['sample'] = new_id
tpm_merge_ensmblcode = tpm_merge_ensmblcode.drop(['index'], axis=1)
first_column = tpm_merge_ensmblcode.pop('sample')
tpm_merge_ensmblcode.insert(0, 'sample', first_column)
tpm_merge_ensmblcode['case'] = tpm_merge_ensmblcode['sample'].str.split('-')
clinical = clinical.rename(columns={'case_submitter_id': 'case'})
clinical = clinical.drop_duplicates(subset='case', keep="last")
length_code=clinical['case'][1].split('-')
tpm_merge_ensmblcode['case'] = tpm_merge_ensmblcode['case'].apply(lambda x: '-'.join(x[:len(length_code)]))

print("Output ready!")
# Generate file paths
tpm_clinical_file_path = os.path.join(folders + name_in_output, name_in_output+'_tpm_clinical.csv')
raw_counts_file_path = os.path.join(folders + name_in_output, name_in_output+'_raw-counts.csv')

print("Saving Clinical data file...")
#merge and save the TMP file with the clinical data
final_tpm_file = pd.merge(tpm_merge_ensmblcode, clinical, on="case")
final_tpm_file.to_csv(tpm_clinical_file_path, index=False)

ensembl_mage= {'ENSG00000198681': 'MAGEA1', 'ENSG00000268606': 'MAGEA2', 'ENSG00000183305': 'MAGEA2B',
				   'ENSG00000221867': 'MAGEA3', 'ENSG00000147381': 'MAGEA4', 'ENSG00000242520': 'MAGEA5P',
				   'ENSG00000197172': 'MAGEA6', 'ENSG00000224732': 'MAGEA7P', 'ENSG00000156009': 'MAGEA8',
				   'ENSG00000123584': 'MAGEA9', 'ENSG00000267978': 'MAGEA9B', 'ENSG00000124260': 'MAGEA10',
				   'ENSG00000185247': 'MAGEA11',
				   'ENSG00000213401': 'MAGEA12', 'ENSG00000214107': 'MAGEB1', 'ENSG00000099399': 'MAGEB2',
				   'ENSG00000198798': 'MAGEB3', 'ENSG00000120289': 'MAGEB4', 'ENSG00000188408': 'MAGEB5',
				   'ENSG00000176746': 'MAGEB6',
				   'ENSG00000177689': 'MAGEB10', 'ENSG00000189023': 'MAGEB16', 'ENSG00000182798': 'MAGEB17',
				   'ENSG00000176774': 'MAGEB18', 'ENSG00000155495': 'MAGEC1', 'ENSG00000046774': 'MAGEC2',
				   'ENSG00000165509': 'MAGEC3',
				   'ENSG00000179222': 'MAGED1', 'ENSG00000102316': 'MAGED2', 'ENSG00000154545': 'MAGED4',
				   'ENSG00000187243': 'MAGED4B', 'ENSG00000198934': 'MAGEE1', 'ENSG00000186675': 'MAGEE2',
				   'ENSG00000177383': 'MAGEF1',
				   'ENSG00000187601': 'MAGEH1', 'ENSG00000254585': 'MAGEL2','ENSG00000268916':'CSAG3',
				   'ENSG00000198930':'CSAG1','ENSG00000268902': 'CSAG2','ENSG00000242599':'CSAG4','ENSG00000130726':'trim28'}
# Rename the columns in the DataFrame based on the dictionary
# Load the CSV file
data = pd.read_csv(tpm_clinical_file_path)
# Rename the columns based on the dictionary
data_renamed = data.rename(columns=ensembl_mage)
# Overwrite the original file with the modified DataFrame
data_renamed.to_csv(tpm_clinical_file_path, index=False)
print("tpm file has been successfully updated with the renamed columns.")
print("Saving raw counts data file...")
#Save the raw counts file
raw_count_table=raw_count_table.set_index(['gene_id'])
raw_count_table=raw_count_table.T
raw_count_table=raw_count_table.reset_index()
raw_count_table=raw_count_table.rename(columns={'index': 'sample'})
raw_count_table.to_csv(raw_counts_file_path, index=False)

# Load the existing configuration from config.json
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

# Update the config dictionary with the new file paths
config['tpm_clinical_file_path'] = tpm_clinical_file_path
config['raw_counts_file_path'] = raw_counts_file_path

# Write the updated configuration back to config.json
with open(config_path, 'w') as config_file:
    json.dump(config, config_file, indent=4)

print("Updated config.json with new file paths.")
print("Done!")
