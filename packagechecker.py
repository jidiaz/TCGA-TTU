import importlib.util
import sys

# List of all the packages and modules to check
packages_to_check = [
    "IPython.display", # display, HTML
    "os",
    "json",
    "pandas",
    "numpy",
    "matplotlib.pyplot",
    "statsmodels.api",
    "seaborn",
    "sklearn.decomposition", # PCA, TruncatedSVD
    "sklearn.preprocessing", # StandardScaler
    "sklearn.pipeline", # make_pipeline
    "sklearn.cluster", # KMeans, AgglomerativeClustering, DBSCAN
    "sklearn.metrics", # silhouette_score, accuracy_score
    "sklearn.model_selection", # train_test_split
    "sklearn.neighbors", # KNeighborsClassifier
    "sklearn.impute", # SimpleImputer
    "lifelines", # KaplanMeierFitter
    "matplotlib.colors", # LinearSegmentedColormap
    "warnings",
    "io",
    "ipywidgets", # FileUpload, Button
    "scipy.stats", # expon
    "datetime" # datetime
]

# Function to check if a package is installed and importable
def check_package(package_name):
    try:
        spec = importlib.util.find_spec(package_name)
        if spec is None:
            return False
    except ImportError:
        return False
    return True

# Check each package in the list and collect the ones that are not installed
not_installed_packages = [package for package in packages_to_check if not check_package(package)]

# Map package import names to pip package names if they differ
pip_package_names = {
    "IPython.display": "ipython",
    "matplotlib.pyplot": "matplotlib",
    "statsmodels.api": "statsmodels",
    "sklearn.decomposition": "scikit-learn",
    "sklearn.preprocessing": "scikit-learn",
    "sklearn.pipeline": "scikit-learn",
    "sklearn.cluster": "scikit-learn",
    "sklearn.metrics": "scikit-learn",
    "sklearn.model_selection": "scikit-learn",
    "sklearn.neighbors": "scikit-learn",
    "sklearn.impute": "scikit-learn",
    "ipywidgets": "ipywidgets",
    "scipy.stats": "scipy",
    "lifelines": "lifelines"
}

# Generate installation commands for missing packages
installation_commands = []
for package in not_installed_packages:
    pip_name = pip_package_names.get(package, package)
    command = f"pip install {pip_name}"
    installation_commands.append(command)

# Write the installation commands to a text file
with open("install_missing_packages.txt", "w") as file:
    for command in installation_commands:
        file.write(command + "\n")

if installation_commands:
    print("Installation commands for missing packages have been written to install_missing_packages.txt")
else:
    print("All packages are installed.")
