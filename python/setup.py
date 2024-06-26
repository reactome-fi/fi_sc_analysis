from setuptools import setup, find_packages

__version__ = '0.1.2'

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name='scpy4reactome',
    version=__version__,
    description='python service for single cell analysis in Reactome',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/reactome-fi/fi_sc_analysis',
    author='Guanming Wu',
    author_email='wu@ohsu.edu',
    license='Apache',
    packages=find_packages(),
    install_requires=[
        # It will be much nicer to use the exact versions to avoid any compatibility issue.
        'jsonrpclib-pelix>=0.4.1',
        # Pick whatever is the latest version
        # 'scanpy>=1.9.8', # Fix to this version since it is used for development
        'networkx>=2.4',
        # 'scvelo>=0.2.2',
        # Most likely 3.6.0 will work
        'matplotlib>=3.6.0', # See this https://stackoverflow.com/questions/24251102/from-matplotlib-import-ft2font-importerror-dll-load-failed-the-specified-pro
        # Use >= under mac to solve h5d issue
        # 'tables>=3.5.2',# Make sure this version is used. Otherwise, there is a dll cannot find error under windows.
        'python-igraph>=0.8.2',
        # 'scikit-learn>=1.4.0',  
        'leidenalg>=0.8.1',
        # 'umap_learn==0.5.1', # This version works with scanpy 1.9.1 and python 3.10
        # 'magic_impute>=2.0.3',
        # 'pyscenic>=0.12.1',
        'gseapy==0.10.4',
        'joblib==0.16.0',
        # For nichenet-based ligand analysis. Make sure nichenetr is installed!!!
        'rpy2>=3.5.12',
        'cellphonedb>=5.0.0',
        'ktplotspy>=0.2.3,' # Use to plot cellphonedb analysis results
        'py4cytoscape>=1.9.0', # Visualize networks in Cytoscape desktop
        'langchain>=0.1.16', # Some LLM stuff
        'langchain_openai>=0.1.3',
        'python-dotenv>=1.0.1'
    ],
    extras_require={
        'pandas': ['pandas>=1.0.5'],
        'numpy': ['numpy>=1.19.1'],
        'scipy': ['scipy>=1.4.1'],
        'anndata': ['anndata>=0.7.4'],
        'statsmodels': ['statsmodels>=0.11.1'],
    },
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.7',
        'License :: OSI Approved :: Apache Software License',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    platforms=['any'],
    python_requires='>=3.7.0',
)
