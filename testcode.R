library(data.tree)
library(DiagrammeR)
library(reshape2)
library(readxl)
library(igraph)

# make sure your working directory is correct
getwd()
setwd('C:/Users/Jingwei Liu/OneDrive - Auburn University/INSY7130/Class materials/BIRCH clustering/BIRCH-Clustering-R-code')
getwd()


source('BIRCH_Source.R')

newfile <- as.data.frame(read_excel(file.choose()))  # Search for the the file
View(newfile)

# Calculate the number of observations of the data. We use this for sampling.
n <- nrow(newfile)
n 

#For example, if we only work with a subset of obs. we draw a sample of size n without replacement from the ID's of the data (i.e. the number of the row)
set.seed(3) #Set the seed to replicate results
s_size = 300 # sample size
sample.id <- sample(1:n, s_size, replace = FALSE) #sampling
sample.id

#Filter only the selected ID's
samplenewfile <- newfile[sample.id, 1:2] # get samples from dataset
View(samplenewfile) #We can view the sample of the dataset


# do BIRCH clustering
tree <- BIRCHCluster(samplenewfile,BranchingFactor = 3, LeafEntries = 6, Threshold = 2)
tree

# you can check what custom field in this tree
# Center is the center of the cluster
# CF is how many CFs in this tree. Only applied in base node(the top node)
# LN is how many leaf nodes in this tree. Only applied in base node(the top node)
# NLN is how many non-leaf nodes in this tree. Only applied in base node(the top node)
# Root is how many root nodes in this tree. Only applied in base node(the top node)
# N is how many points in cluster
# LS is the linear sum in cluster
# SS is square sum in cluster
# List is the index of the observations in cluster
# Radius is the Radius in CF node. 
tree$fieldsAll

# you can check the fields by using print fuction
print(tree,"Center")
print(tree,"Radius")
print(tree, "N", "LS", "SS" )
print(tree,"CF","LN","NLN","Root")

# you can plot the tree
plot(tree)

# plot dendrogram
plot(as.dendrogram(tree))

# a different way to plot the tree
plot(as.igraph(tree, directed = TRUE, direction = "climb"))

# For more visualisations, please see https://cran.r-project.org/web/packages/data.tree/vignettes/data.tree.html#plotting

############################################################################################################################################
# if the data has a certain class, you can check how the observations locate in each cluster
TargetCluster <- FindNode(tree,"LN5")
df <- GetObs(newfile,TargetCluster)
table(df$Class)


# get the all same level node centers to a dataframe. You can use this dataframe for further clustering method
# level is the level of the node: 
#   1 is Base node
#   2 is Root node
#   3 is Non-leaf node
#   4 is Leaf node
#   5 is CF node
CFcenter <- ToDataFrameByLevel(tree,level = 5)
