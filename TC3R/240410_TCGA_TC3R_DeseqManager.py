import os
import json
import csv

# Assuming script_dir is correctly defined as before, fetching the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
config_path = os.path.join(script_dir, 'config.json')

# Load JSON configuration as default settings
with open(config_path, 'r') as config_file:
    config = json.load(config_file)
    deseq_directory = config['deseq_directory']

# Function to create directories and return the output directory path
def create_directories(base_output_dir, gene, experiment):
    output_dir = os.path.join(base_output_dir, gene, experiment)
    os.makedirs(output_dir, exist_ok=True)
    return output_dir

# Function to process files and populate the directory mappings CSV
def process_files(input_dir, output_dir, csv_file_path):
    directory_mappings = []

    # Listing all CSV files in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith(".csv"):
            parts = file_name.rstrip(".csv").split('_')
            gene = parts[3]
            experiment = parts[4]
            file_type = 'DES' if file_name.endswith("DES.csv") else 'RAW'

            # Constructing file paths
            input_file_path = os.path.join(input_dir, file_name)
            output_sub_dir = create_directories(output_dir, gene, experiment)

            # Find or create the directory mapping for this gene and experiment
            mapping = next((m for m in directory_mappings if m['gene'] == gene and m['experiment'] == experiment), None)
            if not mapping:
                mapping = {'gene': gene, 'experiment': experiment, 'input_des': '', 'input_raw': '', 'output_dir': output_sub_dir}
                directory_mappings.append(mapping)
            
            if file_type == 'DES':
                mapping['input_des'] = input_file_path
            else:
                mapping['input_raw'] = input_file_path

    # Writing the directory mappings to a CSV file
    with open(csv_file_path, 'w', newline='') as file:
        fieldnames = ['gene', 'experiment', 'input_des', 'input_raw', 'output_dir']
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        for mapping in directory_mappings:
            writer.writerow(mapping)

# Define the input and output directories and the path for the CSV file
deseq_input_directory = os.path.join(deseq_directory, '240403_TCGA_GBM_DEA_Input')
deseq_output_directory = os.path.join(deseq_directory, '240403_TCGA_GBM_DEA_Output')
mapping_csv_file_path = os.path.join(deseq_output_directory, 'directory_mapping.csv')

# Ensure the output directory exists
os.makedirs(deseq_output_directory, exist_ok=True)

# Process the files and populate the CSV
process_files(deseq_input_directory, deseq_output_directory, mapping_csv_file_path)
