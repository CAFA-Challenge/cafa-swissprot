.. CAFA Swissprot Utilities documentation master file, created by
   sphinx-quickstart on Fri Jan 29 14:23:22 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the CAFA Swissprot Utility documentation!
====================================================

.. danger:: This is rough-draft, preliminary, in-progress documentation for the CAFA Swissprot Utilities package. It is guaranteed to change.

Install the package
-------------------

.. code-block:: shell

   python -m pip install -e "git+https://github.com/CAFA-Challenge/cafa-swissprot#egg=cafa-swissprot"

Download one or more SwissProt files. The most recent release is available here:
ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz

The `uniprot.org` *.gz files contain several files, only the `...sprot-only...` files are relevant here.

Identify low-annotation candidate proteins in a SwissProt file
--------------------------------------------------------------

   1. create a yaml-formatted configuration file. See ``yaml/filter_sprot_species_example.yml`` in the git repo for details
      regarding the contents of the yaml file.

      Here's a minimal example::

        sprot_file: ./uniprot_sprot.dat
        output_directory: ./filtered_output_example
        taxonomies:
           - 9606 # Human
           - 4577 # Zea mays
        ontologies:
           - CCO
           - BPO
           - MFO
        allowed_evidence_codes: [IEA, NR, ND, IC, NAS, TAS, ISS, ISO, ISA, ISM, IGC, IBA, IBD, IKR, IRD, RCA]


   2. run ``CAFA_generate_no_exp_files <YOUR YAML FILEPATH>``
   3. results are written to ``output_directory`` specified in your yaml file

      Output:
        - ``sp_species.[taxon_id].tfa`` FASTA file containing all the proteins in [taxon_id] (replaced by an actual number, e.g. sp_species.9606.tfa for human.)
        - ``sp_species.[taxon_id].MFO.noexp.tfa`` Assuming MFO is specified in config, FASTA file for proteins that are not experimentally annotated in the MFO namespace.
        - ``sp_species.[taxon_id].BPO.noexp.tfa`` See above.
        - ``sp_species.[taxon_id].CCO.noexp.tfa`` See above.
        - ``sp_species.[taxon_id].all.noexp.tfa`` FASTA file containing proteins that are not experimentally annotated in *any* GO ontology.


Count proteins
--------------
1. create a yaml-formatted configuration file. See ``yaml/sprot_growth_example.yml`` in the git repo for details.

   .. note::
      Please note that this yaml file is different from the one used above for ``CAFA_generate_no_exp_files``

   Here is an example yaml configuration::

      input_files:
        # the .dat files need to be downloaded and unpacked from
        # ftp://ftp.uniprot.org/pub/databases/uniprot/
        - ./uniprot/2019_01/uniprot_sprot.dat
        - ./uniprot/2019_06/uniprot_sprot.dat
        - ./uniprot/2020_01/uniprot_sprot.dat
        - ./uniprot/2020_05/uniprot_sprot.dat

      # Taxonomies IDs can be found at
      # https://www.uniprot.org/docs/speclist.txt
      taxonomies:
        - 9606    # Human
        - 4577    # Zea mays
        - 10090   # mouse
        - 10116   # rat
        - 3702    # arabidopsis
        - 7227    # drosophila
        - 559292  # yeast
        - 7955    # zebrafish
        - 44689   # dicty
        - 6239    # worm
        - 83333   # ecoli

      ontologies:
        - CCO
        - BPO
        - MFO

      evidence_codes:
        - EXP
        - IDA
        - IPI
        - IMP
        - IGI
        - IEP


2. run ``CAFA_experimental_growth <YOUR YAML FILEPATH>``
3. results are written to stdout::

     $ CAFA_experimental_growth sprot_growth_example.yml
     filename	taxon	CCO	BPO	MFO
     uniprot/2019_01/uniprot_sprot.dat	9606	9905	8786	7539
     uniprot/2019_06/uniprot_sprot.dat	9606	9996	8960	7533
     uniprot/2020_01/uniprot_sprot.dat	9606	10153	9091	7638
     uniprot/2020_05/uniprot_sprot.dat	9606	10294	9282	8175
     uniprot/2019_01/uniprot_sprot.dat	4577	65      94      80
     uniprot/2019_06/uniprot_sprot.dat	4577	73      111     95
     uniprot/2020_01/uniprot_sprot.dat	4577	73      115     101
     uniprot/2020_05/uniprot_sprot.dat	4577	73      116     102


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   propagation_map
   make_blast_annotation
   test

