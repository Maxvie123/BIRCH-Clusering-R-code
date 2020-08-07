# BIRCH-Clustering-R-code
This is a r package for BIRCH clustering 


## Introduction
This package perform a BIRCH clustering and returns a data.tree structure. 
For more about data.tree, please refer to https://cran.r-project.org/web/packages/data.tree/vignettes/data.tree.html


## About the function
You need to provide 4 inputs to the BIRCH clustering function:
  1. data which is a dataframe that you want to do clustering. 
  2. BranchingFactor which is the maximum children allowed for a non-leaf node.
  3. LeafEntries which is the maximum entries(CFs) allowed for a leaf node.
  4. Threshold which is an upper limit to the radius of a CF.
  
This BIRCH function doesn't contain normalization, if needed, please normalize the data before you use this package.
This BIRCH function is order sensitive which means if you have the same order in your dataframe, you will get the same result.

## About the returned data.tree
There are 10 custom fields in the returned data.tree:
  1. Center is the center of the node (cluster)
  2. CF is how many CFs in this tree. Only applied in base node(the top node)
  3. LN is how many leaf nodes in this tree. Only applied in base node(the top node)
  4. NLN is how many non-leaf nodes in this tree. Only applied in base node(the top node)
  5. Root is how many root nodes in this tree. Only applied in base node(the top node)
  6. N is how many points in cluster
  7. LS is the linear sum in cluster
  8. SS is square sum in cluster
  9. List is the index of the observations in cluster
 10. Radius is the Radius in CF node. 
 
## How to use this package
Please refer to the testcode.R
