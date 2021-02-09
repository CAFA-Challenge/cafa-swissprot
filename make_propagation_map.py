import pandas as pd
from goatools.obo_parser import GODag

def get_cco_propagation_map():
    ''' stub of function get generating a ontology annotation propagation map

    returns a pandas.DataFrame
    '''
    obo_filepath = "/home/scott/Documents/MATLAB/CAFA2/ontology/CAFA3/go_cafa3.obo"
    optional_attrs = ["relationship", "replaced_by", "consider"]
    dag = GODag(
        obo_filepath, optional_attrs=optional_attrs, load_obsolete=False, prt=None
    )

    namespace = 'CCO'
    namespace_long = "cellular_component"

    for id, obj in dag.items():
        if id != obj.item_id:
            print("MISMTACH:")
            print(id)
            print(obj.item_id)
            print(obj.alt_ids)
['a', 'b', 'c']

    propagation_map = {term_object.item_id: [*list(term_object.get_all_upper()), term_object.item_id] for term_id, term_object in dag.items() if term_object.namespace == namespace_long}
    print(len(propagation_map))
    #print(propagation_map)

    df_keys = sorted(list(propagation_map.keys()))

    propagation_map_df = pd.DataFrame(data=0, index=df_keys, columns=df_keys)

    for key, val in propagation_map.items():
        propagation_map_df.loc[key, val] = 1

    return propagation_map_df

if __name__ == "__main__":
    prop_map = get_cco_propagation_map()

    print(prop_map.head())
    print("")
    print(prop_map.index)

    tmp = prop_map.loc["GO:0000112",:]
    print(tmp[tmp == 1].index)

    #prop_map.to_pickle("./cco_propagation_map.pkl")