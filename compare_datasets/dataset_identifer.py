#!/usr/bin/env python3

"""
classifier_dataset_bycount
classify the dataset based on the count (if != 0) in column
:param id1: int 1
:param id2: int 2
:param id1name: output if col1 is not 0 
:param id2name: output if value 2 exists (i.e not NA)
:returns <"Both"> <id1name> <id2name> <"NaN">
"""
def classifier_dataset_bycount(col1,col2,id1name,id2name):
    if col1 != 0 and col2 != 0: return "Both"
    elif col1 == 0 and col2 != 0: return id2name
    elif col1 != 0 and col2 == 0: return id1name
    else: return "NaN"


"""
classify the dataset based on the presence or absence of terms (isoform name) in column
:param id1: value 1
:param id2: value 2
:param id1name: output if value 1 exists (i.e not NA)
:param id2name: output if value 2 exists (i.e not NA)
:returns <"Both"> <id1name> <id2name> <"NaN">
"""
def classifier_dataset(id1, id2, id1name, id2name):
    if not pd.isna(id1) and not pd.isna(id2): return("Both")
    elif pd.isna(id1) and not pd.isna(id2): return(id2name)
    elif not pd.isna(id1) and pd.isna(id2): return(id1name)
    else: return("NaN")
    
    
"""
create a column with the amalgamated isoform id in the abundance file for downstream purposes
currently file has column: isoseq_isoform, ont_isoform with NAs if unique to dataset
apply function after classifier_dataset()
:param dataset: <"Both"> <id1name> <id2name> <"NaN">
:param isoseq_isoform: output if dataset is "Iso-Seq"
:param ont_isoform: output if dataset is "ONT"
:param matched_id: output if dataset is "Both"
:returns <isoseq_isoform"> <ont_isoform> <matched_id> <"NaN">
"""
def unionise_id(dataset, isoseq_isoform, ont_isoform, matched_id):
    if dataset == "Both": return(matched_id)
    elif dataset == "ONT": return(ont_isoform)
    elif dataset == "Iso-Seq": return(isoseq_isoform)
    else: return("NaN")
