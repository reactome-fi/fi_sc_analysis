from setuptools import setup, find_packages

__version__ = '0.0.1'

with open("README.md", "r") as f:
    long_description = f.read()

setup(
    name='name of package',
    version=__version__,
    description='python service for single cell analysis in reactome',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/reactome-fi/fi_sc_analysis',
    author='Guanming Wu',
    author_email='wu@ohsu.edu',
    # which license ?
    license='Apache',
    packages=find_packages(),
    install_requires=[
        'scanpy',
        'anndata',
        'networkx',
        'scvelo',

    ],
    extras_require={
        # fill in versions
        'pandas': ['pandas=='],
        'numpy': ['numpy=='],
        'scipy': ['scipy=='],
    },
    classifiers=[
        'Programming Language :: Python',
        'License :: OSI Approved :: Apache Software License',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: Apache Software License',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    platforms=['any'],
    # which python version
    python_requires='>=3.8',
)