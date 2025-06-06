from setuptools import setup, find_packages

setup(
    name='distanX',
    version='0.2.0',
    packages=find_packages(),
    install_requires=['opencv-python', 'anndata', 'pandas', 'numpy', 'joblib'],
    author='kusurin',
    description='A Python Package for Getting ROIs and Calculating Group Distances in Spatial Transcriptomics Data',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/kusurin/distanX',
    license='AGPL-3.0',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)