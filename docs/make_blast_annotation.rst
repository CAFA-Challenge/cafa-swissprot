#################################
Building CAFA Baseline Predictors
#################################

Naive Predictor
+++++++++++++++

.. code-block:: python

    annotation_filepath = "./input_annotation_for_blast_predictor.txt"
    obo_filepath = "./go_cafa.obo"
    namespace = "cellular_component"

    naive_predictor = generate_naive_predictor(
        annotation_filepath=annotation_filepath,
        obo_filepath=obo_filepath,
        namespace_verbose=namespace,
    )

`generate_naive_predictor()` returns a Pandas Series with GO terms as indices and
floating point numbers as values. Because the naive predictor is based solely on
GO term frequency, there is no protein data in the predictor pandas.Series. Instead,
*every* protein would have a predicted annotation at the given confidence level.

Here is a truncated output example for `generate_naive_predictor()`:

===========   ==========
GO Term       Confidence
===========   ==========
GO\:0000015   0.000000
GO\:0000109   0.000541
GO\:0000110   0.000000
GO\:0000111   0.000000
GO\:0000112   0.000000
GO\:0000113   0.000541
GO\:0000118   0.000000
GO\:0000120   0.000000
GO\:0000123   0.006490
GO\:0000124   0.001082
GO\:0000125   0.000000
GO\:0000126   0.000000
GO\:0000137   0.000000
GO\:0000138   0.000000
GO\:0000139   0.002163
GO\:0000142   0.000541
...           ...
===========   ==========


Blast Predictor
+++++++++++++++

Generating a CAFA baseline Blast annotation requires a few things:

1. An OBO file
2. A Swissprot file
3. A file containing the output of a `blastp` query
4. a file containing a pickled Pandas DataFrame representing the GO ontology hierarchy

Here's a simple example:

.. code-block:: python

   OBO_FILEPATH = './go_cafa.obo'
   SPROT_FILEPATH = './uniprot_sprot.dat'
   BLASTP_SIMILARITY_FILEPATH = './blastp.out'
   PROPAGATION_MAP_FILEPATH = './cco_propagation_map.pkl'
   PROTEINS_OF_INTEREST = {'1A1L2_HUMAN', 'ZZEF1_HUMAN'}
   ONTOLOGIES = {'CCO',}

   blast_generator = BlastPredictorGenerator(
        blastp_results_filepath=BLASTP_SIMILARITY_FILEPATH,
        swissprot_filepath=SPROT_FILEPATH,
        propagation_map_filepath=PROPAGATION_MAP_FILEPATH,
        obo_filepath=OBO_FILEPATH,
   )

   blast_predictions = blast_generator.get_predictions(
      query_protein_ids=PROTEINS_OF_INTEREST
      ontology_namespaces=ONTOLOGIES
   )



`get_predictions()` will return a pandas.DataFrame with proteins as indices (rows)
and GO terms as columns. Values in the dataframe will be floating point numbers from 0 to 1
representing the confidence for the given annotation.

Here's a truncated example of the resulting DataFrame:

===========  ===========  ===========  ===========  ===========    ===
Protein      GO\:0000015  GO\:0000109  GO\:0000111  GO\:0000124    ...
===========  ===========  ===========  ===========  ===========    ===
1A1L2_HUMAN           0           0           0           0        ...
ZZEF1_HUMAN           0           0         0.4         0.4        ...
===========  ===========  ===========  ===========  ===========    ===
