#!/usr/bin/env python
import os
import pandas as pd
import datetime
from functools import reduce
import argparse


#python3 240827_TCGA_datamerge_v5.py -p /Users/juansola/Documents/TC3R/research/TCGA/240919_TCGA-LIHC/240808_TCGA_LIHC_Tumor/ -r 240808_TCGA_LIHC_sample_sheet.tsv -d clinical.tsv -o /Users/juansola/Desktop/prueba/

parser = argparse.ArgumentParser(description='Concatenates TCGA transcriptome raw counts files')
parser.add_argument("-d", required=True, type=argparse.FileType('r'),help='clinical data')
parser.add_argument("-r", required=True, type=argparse.FileType('r'),help='sample sheet from TCGA')
parser.add_argument("-p", required=True,action='store',help='path to the raw count files' )
parser.add_argument("-n", action='store_true', help='A boolean flag to point if the file need to be rename, changing the tcga code for the sample anme')
parser.add_argument("-o", required=True,action='store',help='path to save files' )


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
					files_to_merge.append(z)
				else:
					print("new category"+" "+str(z))
			


else:
	print("No name change")
	for filename in os.listdir(path):
		filename = filename.strip('\n')
		files_to_merge.append(filename)


files_to_merge.remove('.DS_Store')

#Date used to name output files
current_datetime = datetime.datetime.now()
date_format = "%Y%m%d"
formatted_date = current_datetime.strftime(date_format)
date_in_name = formatted_date[2:]



# The below function is used to select just the tpm information for mage genes
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
		'ENSG00000198930','ENSG00000268902','ENSG00000242599','ENSG00000130726','ENSMUSG00000079349','ENSG00000067445',
		'ENSG00000185115','ENSG00000182636','ENSG00000162344','ENSG00000118972','ENSG00000105550','ENSG00000133116','ENSG00000134962','ENSG00000188501']








	list = dataframe.gene_id.tolist()
	for x in list:
		z = x.split(".")[0]
		dataframe.gene_id = dataframe.gene_id.replace(x, z)

	dataframe = dataframe.loc[dataframe['gene_id'].isin(mage_ensmbl)]


	return(dataframe)

# Dcitionary with the relation of ensmelb code and mage gene name
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
				   'ENSG00000198930':'CSAG1','ENSG00000268902': 'CSAG2','ENSG00000242599':'CSAG4','ENSG00000130726':'TRIM28','ENSMUSG00000079349':'MAGEA5',
                   'ENSG00000067445':'TRO','ENSG00000185115':'NSMCE3','ENSG00000182636':'NDN','ENSG00000162344':'FGF19','ENSG00000118972':'FGF23','ENSG00000105550':'FGF21',
				   'ENSG00000133116':'KL','ENSG00000134962':'KLB','ENSG00000188501':'LCTL'}

#If the option to split the data set is selected, the script will concat all files into one single dataframe, labeling the samples based on the threshold.
#The samples will be labeled using the threshold: high (above the threshold), middle (between the threshold and zero), No_expression (TPM == 0)




#extract the name of the project
poject_in_name=clinical['project_id'][0]
orverview_folder = date_in_name + "_" + poject_in_name + '_Overview'

tpm={}
high={}

for x in files_to_merge:
	m=x.strip('\n')
	k=pd.read_csv(path+m,header=None )
	headers=['gene_id','gene_name','gene_type','unstranded','stranded_first','stranded_second','tpm_unstranded','fpkm_unstranded','fpkm_uq_unstranded']
	k=k.iloc[6:]
	k[headers] = k[0].str.split('\t', expand=True)
	k.drop(columns=k.columns[0], axis=1, inplace=True)
	m = m.strip(' ')
	k=k[['gene_id','unstranded', 'tpm_unstranded']]
	k = k.rename(columns={'unstranded': m})
	t=k[['gene_id', 'tpm_unstranded']]
	t=t.rename(columns={'tpm_unstranded': m})
	tpm[m]=t
	high[m]=k[['gene_id', m]]


if os.path.exists(folders + orverview_folder) and os.path.isdir(folders + orverview_folder):
	print("overview file was created")

else:
	os.makedirs(os.path.join(folders + orverview_folder))
	tpm_merge = reduce(lambda left, right: pd.merge(left, right, on=['gene_id'], how='left'), tpm.values())
	raw_merge= reduce(lambda left, right: pd.merge(left, right, on=['gene_id'], how='left'), high.values())

	tpm_merge_ensmblcode = tpmoverview(tpm_merge)


	tpm_merge_ensmblcode = tpm_merge_ensmblcode.set_index(['gene_id'])
	tpm_merge_ensmblcode = tpm_merge_ensmblcode.T
	tpm_merge_ensmblcode = tpm_merge_ensmblcode.reset_index()
	codes = tpm_merge_ensmblcode['index']

	#generate a list with the sample codes we will use to merge the clinical data
	new_id = []
	for x in codes:
		z = x.split('-normal')
		y = z[0].split('-')
		change = y[0] + '-' + y[1] + '-' + y[2]
		new_id.append(change)

	#Merge the tpm values with the clinical information
	tpm_merge_ensmblcode=tpm_merge_ensmblcode.rename(columns={'index':'samples'})
	tpm_merge_ensmblcode['new_id'] = new_id
	tpm_merge_ensmblcode['case'] = tpm_merge_ensmblcode['new_id'].str.split('-')
	clinical = clinical.rename(columns={'case_submitter_id': 'case'})
	clinical = clinical.drop_duplicates(subset='case', keep="last")
	length_code=clinical['case'][1].split('-')
	tpm_merge_ensmblcode['case'] = tpm_merge_ensmblcode['case'].apply(lambda x: '-'.join(x[:len(length_code)]))
	final_tpm_file = pd.merge(tpm_merge_ensmblcode, clinical, on="case")
	name_tpm= final_tpm_file['project_id'][0]
	#label normal tissue
	final_tpm_file = final_tpm_file.drop(['new_id'], axis=1)
	final_tpm_file=final_tpm_file.rename(columns=ensembl_mage)
	#save the files
	final_tpm_file.to_csv(folders + orverview_folder+ '/'+name_tpm+'mage_tpm_plus_clinical.csv',index=False)
	tpm_merge.to_csv(folders + orverview_folder+ '/'+name_tpm+'_all_tpm.csv',index=False)
	raw_merge.to_csv(folders + orverview_folder+ '/'+name_tpm+'_all_raw.csv',index=False)






	