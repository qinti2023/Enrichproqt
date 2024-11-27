from setuptools import setup, find_packages

setup(
    name='Enrichproqt',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'statsmodels'
    ],
    description='Pathway Enrichment Analysis for Protein Lists',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='qinti',
    author_email='qinti@zju.edu.cn',
    url='https://github.com/qinti2023/Enrichproqt', 
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    include_package_data=True,
    package_data={
        'Enrichproqt': ['data/*.pkl'],  
    },
)
