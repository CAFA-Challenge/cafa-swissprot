#################################
Propagation Map Matrices
#################################

Many operations in the codebase rely on GO namespace data in the form of a
Pandas DataFrame. These dataframes contain a binary representation of the hierachy of the ontology and
facilitate term propagation.

Here is example code for generating and saving a propagation matrix:

.. code-block:: python

   obo_filepath = "./cafa.obo"
   namespace = "cellular_component"
   cco_propagation_matrix = get_propagation_map(obo_filepath=obofilepath, namespace=namespace)
   cco_propagation_matrix.to_pickle("./cco_propagation_matrix.pkl")

The propagation matrix DataFrame will look something like this:

===========  =========== =========== =========== =========== =========== ===========
GO Term      GO\:0000015 GO\:0000109 GO\:0000110 GO\:0000111 GO\:0000112 GO\:0000113
===========  =========== =========== =========== =========== =========== ===========
GO\:0000015           1           0           0           0           0  ...
GO\:0000109           0           1           0           0           0  ...
GO\:0000110           0           1           1           0           0  ...
GO\:0000111           0           1           0           1           0  ...
GO\:0000112           0           1           0           0           1  ...
...          ...         ...         ...         ...         ...         ...
===========  =========== =========== =========== =========== =========== ===========
