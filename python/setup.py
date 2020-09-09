from setuptools import setup, find_packages

# https://semver.org/
__version__ = '0.1.0'

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
        'numpy==1.17.0', # tensorflow version dependency
        'scipy==1.4.1',  # tensorflow version dependency
        'tensorflow>=2.3.0',
        'keras>=2.4.3',
        'jsonrpclib-pelix>=0.4.1',
        'scanpy>=1.5.1',
        'networkx>=2.4',
        'scvelo>=0.2.2',
        'python-igraph>=0.8.2',
        'scikit-learn>=0.23.1',
        'leidenalg>=0.8.1',
        'umap_learn>=0.4.6',
        'umap>=0.1.1',
    ],
    extras_require={
        'pandas': ['pandas>=1.0.5'],
        # 'numpy': ['numpy>=1.19.1'],
        # 'scipy': ['scipy>=1.4.1'],
        'anndata': ['anndata>=0.7.4'],
        'statsmodels': ['statsmodels>=0.11.1'],
    },
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3.6',
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
