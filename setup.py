from setuptools import setup, find_packages
version = {}
with open("bactipipe/__version__.py") as f:
    exec(f.read(), version)

setup(
    name="bactipipe",
    version=version["__version__"],
    packages=find_packages(),  # Finds all packages under bactipipe/
    description="Bacterial WGS data analysis pipeline",
    author="Maurice Byukusenge",
    author_email="bmaurice@psu.edu",
    install_requires=[
        "boto3>=1.40.6",
        "tqdm>=4.67",
        "biopython>=1.85",
        "tabulate>=0.9",
    ],
    entry_points={
        "console_scripts": [
            "bactipipe=bactipipe.cli:main"
        ]
    },
    include_package_data=True,
    package_data={
        "bactipipe": ["data/lambda.fasta", "data/pathogenic_bacteria.txt", "data/all_vf_category_map.tsv", "data/dx_vf_category_map.tsv", "data/amr_mutation_map.tsv","data/amr_acquired_map.tsv"],
    },
    python_requires=">=3.7",
    url="https://github.com/mauricebyuka/bactiPipe",  # optional
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
