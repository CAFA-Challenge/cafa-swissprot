from setuptools import setup, find_packages

setup(
    name="cafa-swissprot",
    #name="CAFA Swissprot Targets",
    version="0.0.1",
    author="Iddo Friedberg, Scott Zarecor",
    author_email=None,
    maintainer="Iddo Friedberg, Scott Zarecor",
    # url ='',
    description="Package for identifying and creating CAFA benchmark datasets.",
    # long_description = long_description,
    # long_description_content_type ="text/markdown",
    license="MIT",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "CAFA_experimental_growth = cli:cli_experimental_growth",
            "CAFA_print_annotation_counts = cli:cli_print_annotation_counts",
            "CAFA_generate_no_exp_files = cli:cli_generate_no_exp_files",
            "CAFA_generate_protein_fasta = cli:cli_generate_protein_fasta",
        ]
    },
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    keywords="CAFA",
    install_requires=[
        "biopython",
        "PyYAML",
        "click",
    ],
    zip_safe=False,
)


git filter-branch --force --index-filter "git rm --cached --ignore-unmatch ./tests/test_data/uniprot/*" --prune-empty --tag-name-filter cat -- --all