from setuptools import setup, find_packages

setup(
    name="CAFA Swissprot Targets",
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
            "CAFA_experimental_growth = cli:experimental_growth",
            "CAFA_print_annotation_counts = cli:print_annotation_counts",
            "CAFA_generate_no_exp_files = cli:generate_no_exp_files",
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
    ],
    zip_safe=False,
)
