import os
import json
import pandas as pd
from jinja2 import Template

# Dynamically build the path to the config.json file
script_dir = os.path.dirname(os.path.abspath(__file__))  # Gets the directory where the script is located
config_path = os.path.join(script_dir, 'config.json')  # Builds the path to config.json in the same directory
# Load JSON configuration as default settings
with open(config_path, 'r') as config_file:
    config = json.load(config_file)
    overview_directory = config['overview_directory']

# Directories within the overview directory
data_directory = os.path.join(overview_directory, '240325_TCGA_GBM_ToData')
plots_directory = os.path.join(overview_directory, '240325_TCGA_GBM_ToPlots')

# Content mapping: Specifies where each plot/table belongs
content_mapping = {
    '2.1 Distribution of TPMs Expression': {
        'type': 'table',
        'filename': 'DescriptiveStatistics.csv'
    },
    '2.2 plot_gene_histograms_type1(Type1)': {
        'type': 'plot',
        'filename': 'type1_TPM_Expression_histogram.png'
    },
    '2.3 plot_gene_histograms_type2(Type2)': {
        'type': 'plot',
        'filename': 'type2_TPM_Expression_histogram.png'
    },
    '2.4 Composition of Positive and Negative Tumors by MAGEs': {
        'type': 'plot',
        'filename': 'MAGE_PosNeg_Composition.png'
    },
    '3.1 Correlation Graph': {
        'type': 'plot',
        'filename': 'MAGE_PearsonCorrelation_heatmap.png'
    },
    '3.2 Expression class graph': {
        'type': 'plot',
        'filename': 'Cutoff_Generation_MAGEA1-gene_barplot.png'
    },
    '3.3 Heatmap': {
        'type': 'plot',
        'filename': 'SortedExpression_MAGEA1_Heatmap.png'
    },
    '3.4 Kaplan-Meier Graph': {
        'type': 'plot',
        'filename': 'Kaplan-Meier.png'
    }
}

# Process DescriptiveStatistics.csv
desc_stats_path = os.path.join(data_directory, content_mapping['2.1 Distribution of TPMs Expression']['filename'])
if os.path.exists(desc_stats_path):
    descriptive_statistics_df = pd.read_csv(desc_stats_path, index_col=0)
    content_mapping['2.1 Distribution of TPMs Expression']['content'] = descriptive_statistics_df.T.to_html()

# Process plots
for section, info in content_mapping.items():
    if info['type'] == 'plot':
        file_path = os.path.join(plots_directory, info['filename'])
        if os.path.exists(file_path):
            info['content'] = file_path

# HTML Template
html_template = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>TCGA Data Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; }
        .container { max-width: 800px; margin: auto; }
        img { width: 100%; height: auto; }
        table { width: 100%; border-collapse: collapse; }
        td, th { border: 1px solid #dddddd; text-align: left; padding: 8px; }
        th { background-color: #f2f2f2; }
    </style>
</head>
<body>
<div class="container">
    <h1>TCGA Data Analysis Report</h1>
    {% for section, info in content_mapping.items() %}
    <h2>{{ section }}</h2>
    {% if info.type == 'table' %}
    {{ info.content|safe }}
    {% elif info.type == 'plot' %}
    <img src="{{ info.content }}" alt="{{ section }}">
    {% endif %}
    <p>Lorem ipsum placeholder for authors comments on {{ section }}.</p>
    {% endfor %}
</div>
</body>
</html>
"""

# Rendering HTML content
template = Template(html_template)
html_content = template.render(content_mapping=content_mapping)

# Write to HTML file
output_html_path = 'TCGA_Data_Analysis_Report.html'
with open(output_html_path, 'w') as html_file:
    html_file.write(html_content)

print(f"The HTML report has been generated successfully and saved to {output_html_path}.")
