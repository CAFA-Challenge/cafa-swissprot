CAFA SwissProt Targets
This is a part of the CAFA project.

A collection of functionality related to identifying 
proteins with limited GO annotation from SwissProt files.


## Running the code:
### Install the package: 
`python -m pip install -e "git+https://github.com/CAFA-Challenge/cafa-swissprot#egg=cafa-swissprot"`

Download one or more SwissProt files. The most recent release is available here: 
<ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz>

The `uniprot.org` *.gz files contain several files, only the `...sprot-only...` files are relevant here.


### The package exposes a couple of command-line entrypoints:
1. Identify low-annotation candidate proteins in a SwissProt file:

    - create a yaml-formatted configuration file (see `yaml/filter_sprot_species_example.yml`)
    - run `CAFA_generate_no_exp_files <YOUR YAML FILE>`
    - results are written to `output_directory` specified in your yaml file
      
      Output:
        - `sp_species.[taxon_id].tfa` FASTA file containing all the proteins in [taxon_id] (replaced by an actual number, e.g. sp_species.9606.tfa for human.)
        - `sp_species.[taxon_id].MFO.noexp.tfa` Assuming MFO is specified in config, FASTA file for proteins that are not experimentally annotated in the MFO namespace.
        - `sp_species.[taxon_id].BPO.noexp.tfa` See above.
        - `sp_species.[taxon_id].CCO.noexp.tfa` See above.
        - `sp_species.[taxon_id].all.noexp.tfa` FASTA file containing proteins that are not experimentally annotated in *any* GO ontology.

2. Count proteins 
    - create a yaml-formatted configuration file (see `yaml/sprot_growth_example.yml`). Please note that this yaml 
    file is different from the one used above for `CAFA_generate_no_exp_files`
    - run `CAFA_experimental_growth <YOUR YAML FILE>`
    - results are written to stdout:
        <pre>$ CAFA_experimental_growth sprot_growth_example.yml
        filename	taxon	CCO	BPO	MFO
        /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2019_01/uniprot_sprot.dat	9606	9905	8786	7539
        /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2019_06/uniprot_sprot.dat	9606	9996	8960	7533
        /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2020_01/uniprot_sprot.dat	9606	10153	9091	7638
        /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2020_05/uniprot_sprot.dat	9606	10294	9282	8175
        /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2019_01/uniprot_sprot.dat	4577	65	94	80
        /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2019_06/uniprot_sprot.dat	4577	73	111	95
        /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2020_01/uniprot_sprot.dat	4577	73	115	101
        /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2020_05/uniprot_sprot.dat	4577	73	116	102
        </pre>
    
3. Count proteins with annotation (grouped by ontology namespace) in one or more SwissProt files:
    - create a yaml-formatted configuration file (see `yaml/swissprot_counts_example.yml`). Please note that this yaml 
    file is different from the one used above for `CAFA_generate_no_exp_files`
    - run `CAFA_print_annotation_counts <YOUR YAML FILE>`
    - results are written to stdout:
    
        <pre>
        $ CAFA_print_annotation_counts swissprot_counts_example.yml
        ---------------------------------------------------------------------------------------------------------------------------------------
        | COUNTS OF PROTEINS WITH ANNOTATION                                                                                                  |
        ---------------------------------------------------------------------------------------------------------------------------------------
        | FILE                                                                                    | TAXON | CCO    | BPO    | MFO    | GROWTH |
        ---------------------------------------------------------------------------------------------------------------------------------------
        | /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2019_01/uniprot_sprot_truncated.dat | 9606  |  9,905 |  8,786 |  7,539 |        |
        ---------------------------------------------------------------------------------------------------------------------------------------
        | /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2019_06/uniprot_sprot_truncated.dat | 9606  |  9,996 |  8,960 |  7,533 |   +259 |
        ---------------------------------------------------------------------------------------------------------------------------------------
        | /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2020_01/uniprot_sprot_truncated.dat | 9606  | 10,153 |  9,091 |  7,638 |   +393 |
        ---------------------------------------------------------------------------------------------------------------------------------------
        | /home/scott/virtualenvs/cafa_swissprot/code/uniprot/2020_05/uniprot_sprot_truncated.dat | 9606  | 10,294 |  9,282 |  8,175 |   +869 |
        ---------------------------------------------------------------------------------------------------------------------------------------
        </pre>


