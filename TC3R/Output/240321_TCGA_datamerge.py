#!/usr/bin/env python
import os
import pandas as pd
import datetime
from functools import reduce
import argparse

# example how to run it: 

parser = argparse.ArgumentParser(description='Concatenates TCGA transcriptome raw counts files')
parser.add_argument("-p", required=True,action='store',help='path to the TCGA files' )
parser.add_argument("-d", required=True, type=argparse.FileType('r'),help='clinical data')
parser.add_argument("-r", required=True, type=argparse.FileType('r'),help='sample sheet from TCGA')
parser.add_argument("-n", action='store_true', help='A boolean flag to point if the file need to be rename, changing the tcga code for the sample anme')
parser.add_argument("-o", required=True,action='store',help='path to save files' )

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
    print("No name cahnge")


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

for filename in os.listdir(path):
	
	 #Load each file
	filename=filename.strip('\n')
	tcga_file=pd.read_csv(path+filename,header=None)
	
    #remove the rows with information that is not relevant to us - metric about the mapping reads- 
	headers=['gene_id','gene_name','gene_type','unstranded','stranded_first','stranded_second','tpm_unstranded','fpkm_unstranded','fpkm_uq_unstranded']
	tcga_file=tcga_file.iloc[2:]
	tcga_file[headers] = tcga_file[0].str.split('\t', expand=True)
	tcga_file.drop(columns=tcga_file.columns[0], axis=1, inplace=True)
	
    #Select the rawcounts (unstrated column) and tpm values (tpm_unstranded column) and then save the information in a dictionary usign as key the name of the file
	filename = filename.strip(' ')
	tcga_file=tcga_file[4:][['gene_id','unstranded', 'tpm_unstranded']]
	tcga_file= tcga_file.rename(columns={'unstranded': filename})
	tpm=tcga_file[['gene_id', 'tpm_unstranded']]
	tpm=tpm.rename(columns={'tpm_unstranded': filename})
	tpm_clinical[filename]=tpm
	raw_count[filename]=tcga_file[['gene_id', filename]]
	

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

#Select only the genes we are interested on
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

#modifie the sample ID (delete the numer of the replicate) to match with the ID in the clinical file
tpm_merge_ensmblcode['sample'] = new_id
tpm_merge_ensmblcode = tpm_merge_ensmblcode.drop(['index'], axis=1)
first_column = tpm_merge_ensmblcode.pop('sample')
tpm_merge_ensmblcode.insert(0, 'sample', first_column)
tpm_merge_ensmblcode['case'] = tpm_merge_ensmblcode['sample'].str.split('-')
clinical = clinical.rename(columns={'case_submitter_id': 'case'})
clinical = clinical.drop_duplicates(subset='case', keep="last")
length_code=clinical['case'][1].split('-')
tpm_merge_ensmblcode['case'] = tpm_merge_ensmblcode['case'].apply(lambda x: '-'.join(x[:len(length_code)]))

#merge and save the TMP file with the clinical data
final_tpm_file = pd.merge(tpm_merge_ensmblcode, clinical, on="case")
final_tpm_file.to_csv(folders + name_in_output+ '/'+name_in_output+'_tpm_clinical.csv',index=False)

#Save the raw caounts file
raw_count_table=raw_count_table.set_index(['gene_id'])
raw_count_table=raw_count_table.T
raw_count_table=raw_count_table.reset_index()
raw_count_table=raw_count_table.rename(columns={'index': 'sample'})
raw_count_table.to_csv(folders + name_in_output+ '/'+name_in_output+'_raw-counts.csv',index=False)
