#!/usr/bin/env python
import os
import pandas as pd
import datetime
from functools import reduce
import argparse
from collections import Counter
import json

# Dynamically build the path to the config.json file
script_dir = os.path.dirname(os.path.abspath(__file__))  # Gets the directory where the script is located
config_path = os.path.join(script_dir, 'config.json')  # Builds the path to config.json in the same directory
# Load JSON configuration as default settings
with open(config_path, 'r') as config_file:
    config = json.load(config_file)

# example how to run it: 

parser = argparse.ArgumentParser(description='Concatenates TCGA transcriptome raw counts files')
parser.add_argument("-p", default=config['path_to_TCGA_files'], action='store', help='path to the TCGA files')
parser.add_argument("-d", default=config['clinical_data_path'], type=argparse.FileType('r'), help='clinical data')
parser.add_argument("-r", default=config['sample_sheet_path'], type=argparse.FileType('r'), help='sample sheet from TCGA')
parser.add_argument("-n", action='store_true', default=config['rename_files'], help='A boolean flag to point if the file needs to be renamed, changing the TCGA code for the sample name')
parser.add_argument("-o", default=config['output_path'], action='store', help='path to save files')

args = parser.parse_args()
clinical=pd.read_csv(args.d,sep='\t')
sampleinfo=pd.read_csv(args.r,sep='\t')
path=args.p
folders=args.o
name_change=args.n




def rename_duplicates_in_column(df, column_name):
    # Create a counter dictionary to track occurrences
	counts = {}
    # New column to store updated names
	new_names = []
	for name in df[column_name]:
		name= name.split(",")[0]
		if name in counts:
			counts[name] += 1
			new_name = str(name)+"-"+str(counts[name])
		else:
			counts[name] = 1
			new_name = name
		new_names.append(new_name)
    
    # Assign the list of new names to the original column
	df[column_name] = new_names
	return df

# Rename duplicates in the 'Sample ID' column
sampleinfo = rename_duplicates_in_column(sampleinfo, 'Sample ID')


id=[]
files_to_merge=[]
if name_change:
	for filename in os.listdir(path):
		filename = filename.strip('\n')
		id.append(filename)

	for iname in id:
		siname = iname.split('.rna')[0]
		for fname in sampleinfo.index:
			code  = sampleinfo.at[fname, 'File Name']
			sample_code = code.split('.rna_seq')[0]
			if siname == sample_code:
				if "Tumor" in sampleinfo.at[fname, 'Sample Type']:
					z = sampleinfo.at[fname, 'Sample ID']
					os.rename(path+iname, path+z)
					files_to_merge.append(z)
				elif "Normal" in sampleinfo.at[fname, 'Sample Type']:
					z = sampleinfo.at[fname, 'Sample ID']
					z=str(z)+str("-normal")
					os.rename(path+iname, path+z)
				else:
					print("new category"+" "+str(z))
			


else:
	for filename in os.listdir(path):
		filename = filename.strip('\n')
		if "normal" not in filename:
			files_to_merge.append(filename)


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
for filename in files_to_merge:
    # Skip .DS_Store files
    if filename == '.DS_Store':
        skipped_files.append(filename)
        continue

    try:
        # Load each file
        filename = filename.strip('\n')
        tcga_file = pd.read_csv(path + filename, header=None)
        tcga_file = tcga_file.iloc[6:]
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

#Merge files
print("initializing merging process, this may take few minutes...")
#Merge files
raw_count_table = reduce(lambda left, right: pd.merge(left, right, on=['gene_id'], how='left'), raw_count.values())
tpm_count_table = reduce(lambda left, right: pd.merge(left, right, on=['gene_id'], how='left'), tpm_clinical.values())
	

#Function to recover genes of interest 
def tpmoverview (dataframe):
	mage_ensmbl = ['ENSG00000163565', 'ENSG00000175793', 'ENSG00000078018', 'ENSG00000107438', 'ENSG00000138685',
    'ENSG00000163661', 'ENSG00000155324', 'ENSG00000088726', 'ENSG00000137752', 'ENSG00000168542',
    'ENSG00000178031', 'ENSG00000156453', 'ENSG00000182836', 'ENSG00000101916', 'ENSG00000204397',
    'ENSG00000120337', 'ENSG00000138134', 'ENSG00000186832', 'ENSG00000154127', 'ENSG00000107796',
    'ENSG00000165175', 'ENSG00000164904', 'ENSG00000163568', 'ENSG00000139946', 'ENSG00000161671',
    'ENSG00000170458', 'ENSG00000128274', 'ENSG00000008118', 'ENSG00000255221', 'ENSG00000026103',
    'ENSG00000168702', 'ENSG00000150782', 'ENSG00000138271', 'ENSG00000114948', 'ENSG00000106571',
    'ENSG00000138814', 'ENSG00000198835', 'ENSG00000155307', 'ENSG00000163739', 'ENSG00000115884',
    'ENSG00000164035', 'ENSG00000179046', 'ENSG00000169851', 'ENSG00000136237', 'ENSG00000158458',
    'ENSG00000242265', 'ENSG00000276289', 'ENSG00000079308', 'ENSG00000112769', 'ENSG00000133488',
    'ENSG00000103647', 'ENSG00000197594', 'ENSG00000162409', 'ENSG00000156510', 'ENSG00000152766',
    'ENSG00000146678', 'ENSG00000108379', 'ENSG00000177283', 'ENSG00000120251', 'ENSG00000005884',
    'ENSG00000154162', 'ENSG00000143127', 'ENSG00000057294', 'ENSG00000124762', 'ENSG00000111341',
    'ENSG00000198829', 'ENSG00000196878', 'ENSG00000049089', 'ENSG00000145358', 'ENSG00000136546',
    'ENSG00000164692', 'ENSG00000075651', 'ENSG00000117461', 'ENSG00000104267', 'ENSG00000188452',
    'ENSG00000164100', 'ENSG00000163359', 'ENSG00000143578', 'ENSG00000152056', 'ENSG00000130635',
    'ENSG00000074370', 'ENSG00000138449', 'ENSG00000104213', 'ENSG00000164690', 'ENSG00000162552',
    'ENSG00000152661', 'ENSG00000065618', 'ENSG00000177414', 'ENSG00000143375', 'ENSG00000138759',
    'ENSG00000137757', 'ENSG00000157613', 'ENSG00000092969', 'ENSG00000120471', 'ENSG00000137747',
    'ENSG00000128591', 'ENSG00000116711', 'ENSG00000157240', 'ENSG00000162365', 'ENSG00000067445',
    'ENSG00000164651', 'ENSG00000166268', 'ENSG00000166863', 'ENSG00000165168', 'ENSG00000144229',
    'ENSG00000150628', 'ENSG00000182871', 'ENSG00000152315', 'ENSG00000110042', 'ENSG00000134627',
    'ENSG00000197915', 'ENSG00000114251', 'ENSG00000224383', 'ENSG00000249212', 'ENSG00000139329',
    'ENSG00000186469', 'ENSG00000188064', 'ENSG00000115596', 'ENSG00000061492']

	list = dataframe.gene_id.tolist()
	for x in list:
		z = x.split(".")[0]
		dataframe.gene_id = dataframe.gene_id.replace(x, z)

	dataframe = dataframe.loc[dataframe['gene_id'].isin(mage_ensmbl)]


	return(dataframe)

#Select only the genes we are interested on
tpm_merge_ensmblcode = tpmoverview(tpm_count_table)


#
tpm_merge_ensmblcode = tpm_merge_ensmblcode.set_index(['gene_id'])
tpm_merge_ensmblcode = tpm_merge_ensmblcode.T
tpm_merge_ensmblcode = tpm_merge_ensmblcode.reset_index()
codes = tpm_merge_ensmblcode['index']



#modify the sample ID (delete the number of the replicate) to match with the ID in the clinical file
tpm_merge_ensmblcode['sample'] = codes
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

ensembl_mage= {'ENSG00000163565': 'IFI16', 'ENSG00000175793': 'SFN', 'ENSG00000078018': 'MAP2', 'ENSG00000107438': 'PDLIM1', 'ENSG00000138685': 'FGF2',
    'ENSG00000163661': 'PTX3', 'ENSG00000155324': 'GRAMD2B', 'ENSG00000088726': 'TMEM40', 'ENSG00000137752': 'CASP1', 'ENSG00000168542': 'COL3A1',
    'ENSG00000178031': 'ADAMTSL1', 'ENSG00000156453': 'PCDH1', 'ENSG00000182836': 'PLCXD3', 'ENSG00000101916': 'TLR8', 'ENSG00000204397': 'CARD16',
    'ENSG00000120337': 'TNFSF18', 'ENSG00000138134': 'STAMBPL1', 'ENSG00000186832': 'KRT16', 'ENSG00000154127': 'UBASH3B', 'ENSG00000107796': 'ACTA2',
    'ENSG00000165175': 'MID1IP1', 'ENSG00000164904': 'ALDH7A1', 'ENSG00000163568': 'AIM2', 'ENSG00000139946': 'PELI2', 'ENSG00000161671': 'EMC10',
    'ENSG00000170458': 'CD14', 'ENSG00000128274': 'A4GALT', 'ENSG00000008118': 'CAMK1G', 'ENSG00000255221': 'CARD17', 'ENSG00000026103': 'FAS',
    'ENSG00000168702': 'LRP1B', 'ENSG00000150782': 'IL18', 'ENSG00000138271': 'GPR87', 'ENSG00000114948': 'ADAM23', 'ENSG00000106571': 'GLI3',
    'ENSG00000138814': 'PPP3CA', 'ENSG00000198835': 'GJC2', 'ENSG00000155307': 'SAMSN1', 'ENSG00000163739': 'CXCL1', 'ENSG00000115884': 'SDC1',
    'ENSG00000164035': 'EMCN', 'ENSG00000179046': 'TRIML2', 'ENSG00000169851': 'PCDH7', 'ENSG00000136237': 'RAPGEF5', 'ENSG00000158458': 'NRG2',
    'ENSG00000242265': 'PEG10', 'ENSG00000276289': 'KCNE1B', 'ENSG00000079308': 'TNS1', 'ENSG00000112769': 'LAMA4', 'ENSG00000133488': 'SEC14L4',
    'ENSG00000103647': 'CORO2B', 'ENSG00000197594': 'ENPP1', 'ENSG00000162409': 'PRKAA2', 'ENSG00000156510': 'HKDC1', 'ENSG00000152766': 'ANKRD22',
    'ENSG00000146678': 'IGFBP1', 'ENSG00000108379': 'WNT3', 'ENSG00000177283': 'FZD8', 'ENSG00000120251': 'GRIA2', 'ENSG00000005884': 'ITGA3',
    'ENSG00000154162': 'CDH12', 'ENSG00000143127': 'ITGA10', 'ENSG00000057294': 'PKP2', 'ENSG00000124762': 'CDKN1A', 'ENSG00000111341': 'MGP',
    'ENSG00000198829': 'SUCNR1', 'ENSG00000196878': 'LAMB3', 'ENSG00000049089': 'COL9A2', 'ENSG00000145358': 'DDIT4L', 'ENSG00000136546': 'SCN7A',
    'ENSG00000164692': 'COL1A2', 'ENSG00000075651': 'PLD1', 'ENSG00000117461': 'PIK3R3', 'ENSG00000104267': 'CA2', 'ENSG00000188452': 'CERKL',
    'ENSG00000164100': 'NDST3', 'ENSG00000163359': 'COL6A3', 'ENSG00000143578': 'CREB3L4', 'ENSG00000152056': 'AP1S3', 'ENSG00000130635': 'COL5A1',
    'ENSG00000074370': 'ATP2A3', 'ENSG00000138449': 'SLC40A1', 'ENSG00000104213': 'PDGFRL', 'ENSG00000164690': 'SHH', 'ENSG00000162552': 'WNT4',
    'ENSG00000152661': 'GJA1', 'ENSG00000065618': 'COL17A1', 'ENSG00000177414': 'UBE2U', 'ENSG00000143375': 'CGN', 'ENSG00000138759': 'FRAS1',
    'ENSG00000137757': 'CASP5', 'ENSG00000157613': 'CREB3L1', 'ENSG00000092969': 'TGFB2', 'ENSG00000120471': 'TP53AIP1', 'ENSG00000137747': 'TMPRSS13',
    'ENSG00000128591': 'FLNC', 'ENSG00000116711': 'PLA2G4A', 'ENSG00000157240': 'FZD1', 'ENSG00000162365': 'CYP4A22', 'ENSG00000067445': 'TRO',
    'ENSG00000164651': 'SP8', 'ENSG00000166268': 'MYRFL', 'ENSG00000166863': 'TAC3', 'ENSG00000165168': 'CYBB', 'ENSG00000144229': 'THSD7B',
    'ENSG00000150628': 'SPATA4', 'ENSG00000182871': 'COL18A1', 'ENSG00000152315': 'KCNK13', 'ENSG00000110042': 'DTX4', 'ENSG00000134627': 'PIWIL4',
    'ENSG00000197915': 'HRNR', 'ENSG00000114251': 'WNT5A', 'ENSG00000224383': 'PRR29', 'ENSG00000249212': 'ATP1B1P1', 'ENSG00000139329': 'LUM',
    'ENSG00000186469': 'GNG2', 'ENSG00000188064': 'WNT7B', 'ENSG00000115596': 'WNT6', 'ENSG00000061492': 'WNT8A'}
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
#with open(config_path, 'r') as config_file:
#    config = json.load(config_file)

# Update the config dictionary with the new file paths
#config['tpm_clinical_file_path'] = tpm_clinical_file_path
#config['raw_counts_file_path'] = raw_counts_file_path

# Write the updated configuration back to config.json
#with open(config_path, 'w') as config_file:
#    json.dump(config, config_file, indent=4)

#print("Updated config.json with new file paths.")
print("Done!")
