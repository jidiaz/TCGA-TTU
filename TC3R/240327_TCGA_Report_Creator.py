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
data_directory = os.path.join(overview_directory, '240328_TCGA_GBM_Data')
plots_directory = os.path.join(overview_directory, '240328_TCGA_GBM_Plots')

# Content mapping: Specifies where each plot/table belongs
content_mapping = {
    '1.1 Project Overview': {
        'type': 'table',
        'filename': 'ProjectOverview.csv'
    },
    '1.2 Clinical Data Summary': {
        'type': 'table',
        'filename': 'ClinicalDataSummary.csv'
    },       
    '2.1 Distribution of TPMs Expression': {
        'type': 'table',
        'filename': 'DescriptiveStatistics.csv',
        'author_comments': "This table provides an overview of the project."
    },
    '2.2 Type1 Gene Histograms': {
        'type': 'plot',
        'filename': 'type1_TPM_Expression_histogram.png',
        'author_comments': "This plot provides an overview of the project."
    },
    '2.3 Type2 Gene Histograms': {
        'type': 'plot',
        'filename': 'type2_TPM_Expression_histogram.png'
    },
    '2.4 Composition of Positive and Negative Tumors by MAGEs': {
        'type': 'plot',
        'filename': 'MAGE_PosNeg_Composition.png'
    },
    '3.1 Correlation Analysis Summary': {
        'type': 'plot',
        'filename': 'CorrelationAnalysisSummary.csv'
    },
    '3.2 Correlation Graph': {
        'type': 'plot',
        'filename': 'MAGE_PearsonCorrelation_heatmap.png'
    },
    '3.3 Top 10 Correlation Genes': {
        'type': 'table',
        'filename': 'Top_10_Correlation.csv'
    },
    '3.4 Bottom 10 Correlation Genes': {
        'type': 'table',
        'filename': 'Bot_10_Correlation.csv'
    },
    '4.1 Cutoff Generation Summary': {
        'type': 'table',
        'filename': ' CutoffGenerationSummary.csv'
    },
    '4.2 Expression Class Graph': {
        'type': 'plot',
        'filename': "240328_TCGA_GBM_Cutoff/Cutoff_Generation_['MAGEA11']-gene_barplot.png"
    },
    '4.3 Heatmap': {
        'type': 'plot',
        'filename': 'SortedExpression_MAGEA1_Heatmap.png'
    },
    '4.4 Kaplan-Meier Graph': {
        'type': 'plot',
        'filename': 'Kaplan-Meier.png'
    },
    '5.1 Deseq Summary': {
        'type': 'table',
        'filename': 'DeseqSummary.csv'
    }
}

# Process DescriptiveStatistics.csv
desc_stats_path = os.path.join(data_directory, content_mapping['2.1 Distribution of TPMs Expression']['filename'])
if os.path.exists(desc_stats_path):
    descriptive_statistics_df = pd.read_csv(desc_stats_path, index_col=0)
    content_mapping['2.1 Distribution of TPMs Expression']['content'] = descriptive_statistics_df.T.to_html()

# Placeholder text for missing content
placeholder_text = "Content unavailable at this time."

# Adjusted processing logic for tables and plots, including optional author comments
for section, info in content_mapping.items():
    file_path = os.path.join(plots_directory if info['type'] == 'plot' else data_directory, info['filename'])
    author_comments = info.get('author_comments', '')

    # Initialize content with author comments if they exist; different handling for plots and tables
    if author_comments and info['type'] == 'table':
        # For tables, prepend author comments
        content = f"{author_comments}\n"
    elif author_comments and info['type'] == 'plot':
        # For plots, author comments will be appended after processing the file path
        content = ""
    else:
        # No author comments provided
        content = ""

    if os.path.exists(file_path):
        if info['type'] == 'table':
            df = pd.read_csv(file_path, index_col=0)
            info['content'] = f"{content}{df.to_html()}"  # Include author comments above tables if they exist
        elif info['type'] == 'plot':
            # For plots, include the file path and append author comments if they exist
            info['content'] = f"{file_path}\n{content}".strip()  # Ensure we strip in case there are no author comments
    else:
        # Use placeholder text for missing content, without author comments
        info['content'] = placeholder_text

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
