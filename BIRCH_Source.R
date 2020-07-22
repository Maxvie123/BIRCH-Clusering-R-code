
#################################################################################################
# This code perform a BIRCH clustering.  
# To use this code, you need to source this file.
#
# Author: Jingwei Liu
# Version: V1.0
# Last update: 7/23/2020
################################################################################################



# main function
# data is the data you want to do clustering
# BrachingFactor is the maximum children allowed for non-leaf node
# LeafEntries is the maximum entries(CFs) allowed for leaf node
# Threshold is an upper limit to the radius of cluster in CF
BIRCHCluster = function(data, BrachingFactor, LeafEntries,Threshold){
  B <- BrachingFactor
  L <- LeafEntries
  T <- Threshold
  df <- as.data.frame(data)
  ObsNum <- nrow(df)
  tree <- InitialTree(df[1,])
  for (i in 2:ObsNum){
    BuildTree(df[i,],tree,B,L,T)
  }
  return(tree)
}


# create the initial tree
# firstpoint is the row values of one observation
InitialTree = function(firstpoint){
  # Initial the tree, add base node, 1st root, non-leaf, leaf and cf node in the tree
  BirchTree<- Node$new("Tree")
  Root1<- BirchTree$AddChild("Root1")
  NLN1 <- Root1$AddChild("NLN1")
  LN1 <- NLN1$AddChild("LN1")
  CF1 <- LN1$AddChild("CF1")
  
  # Tree node has 9 attributes: The number of points; 
  # The The number of CF nodes; the number of leaf nodes; 
  # the number of non-leaf nodes; the number of root nodes
  # linear sum; square sum; center value ; point index list
  BirchTree$N <-1
  BirchTree$CF <-1
  BirchTree$LN <-1
  BirchTree$NLN <-1
  BirchTree$Root <-1 
  BirchTree$LS <-0
  BirchTree$SS <-0
  BirchTree$Center <- as.numeric(firstpoint)
  BirchTree$List <- c(rownames(firstpoint))
  
  # Root node has 5 attributes: The number of points; linear sum;
  # square sum; center value ; point index list
  Root1$N <-1
  Root1$LS <- 0
  Root1$SS <- 0
  Root1$Center <- as.numeric(firstpoint)
  Root1$List <-c(rownames(firstpoint))
  
  # Non-leaf node has 5 attributes: The number of points; linear sum;
  # square sum; center value; point index list
  NLN1$N <-1
  NLN1$LS <- 0
  NLN1$SS <- 0
  NLN1$Center <-as.numeric(firstpoint)
  NLN1$List <-c(rownames(firstpoint))
  
  # leaf node has 5 attributes: The number of points; linear sum;
  # square sum; center value; point index list
  LN1$N <-1
  LN1$LS <-0
  LN1$SS <- 0
  LN1$Center <-as.numeric(firstpoint)
  LN1$List <-c(rownames(firstpoint))
  
  
  # CF node has 6 attributes: The number of points; center value;
  # point index list; radius; linear sum; square Sum
  CF1$Center <- as.numeric(firstpoint)
  CF1$N <-1
  CF1$LS <-as.numeric(firstpoint)
  CF1$SS <-as.numeric(firstpoint*firstpoint)
  CF1$List <-c(rownames(firstpoint))
  CF1$Radius <- ComputeRadius(CF1$N, CF1$LS, CF1$SS)

  
  return(BirchTree)
}


# build the tree
# datapoint is the row values of one observation.
# tree is the data.tree structure
# BrachingFactor is the maximum children allowed for non-leaf node
# LeafEntries is the maximum entries(CFs) allowed for leaf node
# Threshold is an upper limit to the radius of cluster in CF
BuildTree = function(datapoint, tree, BrachingFactor, LeafEntries,Threshold){

  # find the target CF
  BaseNode <- FindNode(tree, "Tree")
  RootName <- ClosestChild(datapoint,BaseNode)
  RootNode <- FindNode(BaseNode, RootName)
  NLNName <- ClosestChild(datapoint,RootNode)
  NLNNode <- FindNode(RootNode, NLNName)
  LNName <- ClosestChild(datapoint,NLNNode)
  LNNode <- FindNode(NLNNode, LNName)
  CFName <- ClosestChild(datapoint,LNNode)
  CFNode <- FindNode(LNNode, CFName)
  
  # assign datapoint to CF or create new CF
  if (CheckRadius(datapoint, CFNode, Threshold)){
    AssignDataPointToCF(datapoint,CFNode)
  } else{
    CreateNewCF(datapoint, LNNode)
  }
  
  # check if split leaf needed
  if (LNNode$count > LeafEntries){
    SplitNode(LNNode)
  }
  
  # check if split non-leaf needed
  if(NLNNode$count > BrachingFactor ){
    SplitNode(NLNNode)
  }
  
  # check if split root needed
  if (RootNode$count > BrachingFactor){
    SplitNode(RootNode)
  }
  
  # update attributes value
  UpdateAttr(tree)
}



# find the closest child node
# datapoint is the row values of one observation.
# PN is the parent node in the tree. 
ClosestChild = function(datapoint,PN){
  
  # find the centers of all child nodes
  df <- ToDataFrameNetwork(PN, "Center")
  mask <- df$from == PN$name
  centerstring <- df[mask,]$Center
  dimension <- length(datapoint)
  
  # find the closest child node
  flag <- Inf
  for(i in 1:length(centerstring)){
    
    centervector <- vector()
    for (j in 1:dimension){
    value <- as.numeric(strsplit(centerstring[i],",")[[1]][j])
    centervector <- append(centervector, value) 
    }
    euclidean <- dist(rbind(datapoint,centervector))
    if (euclidean[1] < flag){
      flag <- euclidean
      result <- i
    }
  }
  
  # return child node name
  return(PN$children[[result]]$name)
  
}


# Function to compute the CF radius
# N is the number of points
# LS is the linear sum
# SS is the square sum
ComputeRadius = function(N,LS,SS){
  if (N == 1){
    return(0)
  } else{
    X = SS - ((LS*LS)/N)
    R = sqrt(sum(X)/N)
    return(R)
  }
}


# check Threshold constraint
# datapoint is the row values of one observation
# CFNode is the CF node
# Threshold is an upper limit to the radius of cluster in CF
CheckRadius = function(datapoint, CFNode, Threshlod){
  N <-CFNode$N+1
  LS <- CFNode$LS + as.numeric(datapoint)
  SS <- CFNode$SS + as.numeric(datapoint*datapoint)
  if (ComputeRadius(N,LS,SS)<=Threshlod){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# Assign datapoint to the corresponding CF
# datapoint is the row values of one observation you want to add to the corresponding CF
# CFNode is the corresponding CF node
AssignDataPointToCF = function(datapoint, CFNode){
  CFNode$N <- CFNode$N + 1
  CFNode$LS <- CFNode$LS + as.numeric(datapoint)
  CFNode$SS <- CFNode$SS + as.numeric(datapoint*datapoint)
  CFNode$List <- append(CFNode$List, rownames(datapoint))
  CFNode$Radius <- ComputeRadius(CFNode$N,CFNode$LS,CFNode$SS)
  CFNode$Center <- CFNode$LS/CFNode$N
  return(CFNode)
}


# create a new CF node
# datapoint is the row values of one observation you want to add to new CF
# PN is the parent Node of the new CF
CreateNewCF = function(datapoint, PN){
  
  #update root atrribute
  rootNode <- PN$root
  rootNode$CF <- rootNode$CF+1
  
  # Create new CF
  CFName <- paste("CF", PN$root$CF, sep="")
  NewCF <- PN$AddChild(CFName)
  
  # initialize CF attribute
  NewCF$Center <- as.numeric(datapoint)
  NewCF$N <-1
  NewCF$LS <-as.numeric(datapoint)
  NewCF$SS <-as.numeric(datapoint*datapoint)
  NewCF$List <-c(rownames(datapoint))
  NewCF$Radius <- ComputeRadius(NewCF$N, NewCF$LS, NewCF$SS)
  
}

# split Node
# Node is the Node that need to be split
SplitNode = function(Node){
  
  # create new node under parent node
  PN <- Node$parent
  level <- Node$level
  rootNode <- Node$root
  if (level == 4){
    rootNode$LN <- rootNode$LN+1
    prefix = "LN"
    suffix = rootNode$LN
  } else if(level == 3){
    rootNode$NLN <- rootNode$NLN+1
    prefix = "NLN"
    suffix = rootNode$NLN
  } else {
    rootNode$Root <- rootNode$Root+1
    prefix = "Root"
    suffix = rootNode$Root
  }
  NewName <- paste(prefix, suffix, sep = "")
  NodeClone <- Clone(Node)
  NodeClone$name <- NewName
  NewNode <- PN$AddChildNode(NodeClone)
  
  
  # get brach number
  N <- Node$count
  
  # create the center vector of all the child branches
  distvector <- vector()
  for (i in 1:N){
    distvector <- rbind(distvector, Node$children[[i]]$Center)
  }
  
  # create the distance matrix
  distmatrix <- dist(distvector)
  
  # find the largest distance in the matrix and the seeds
  df <- melt(as.matrix(distmatrix),varnames = c("row", "col"))
  target <- df[df$value == max(distmatrix) & df$row>df$col,][1,]
  seed1 <- target$row
  seed2 <- target$col
  
  # form the group for each seed
  seed1vec <- vector()
  seed2vec <- vector()
  seed1vec <- append(seed1vec,Node$children[[seed1]]$name)
  seed2vec <- append(seed2vec,Node$children[[seed2]]$name)
  for (i in 1:N){
    if (i == seed1 | i == seed2 ){
    
    } else{
        dist1 <- df[df$row == seed1&df$col == i,]$value
        dist2 <- df[df$row == seed2&df$col == i,]$value
        if (dist1 >= dist2){
          seed2vec <- append(seed2vec, Node$children[[i]]$name)
        } else{
          seed1vec <- append(seed1vec, Node$children[[i]]$name)
        }
    }

  }
  
  # assign groups to new nodes
  Prune(Node,function(x) (x$name %in% seed1vec)|(x$level>level+1))
  Prune(NewNode,function(x) (x$name %in% seed2vec)|(x$level>level+1))
  
  
}


# update attribute values in the tree
# Node is the root Node of a tree
UpdateAttr = function(Node){
  
  ChildNum <- Node$count
  
  if (Node$level != 5){
    N <- 0
    LS <- 0
    SS <- 0
    List <- vector()
    for (i in 1:ChildNum){
      Child <- Node$children[[i]]
      ChildAttr <- UpdateAttr(Child)
      N <- N + ChildAttr[[1]]
      LS <- LS + ChildAttr[[2]]
      SS <- SS + ChildAttr[[3]]
      List <- append(List, ChildAttr[[4]])
    }
    Node$N <- N
    Node$LS <- LS
    Node$SS <- SS
    Node$Center <- LS/N
    Node$List <- List
    return(list(N,LS,SS,List))

  }else{
    return(list(Node$N, Node$LS, Node$SS, Node$List))
  }

  
  
}


# get the points contained in the node(cluster)
# df is the ORIGINAL dataframe without sampling
# Node is the node contains the points indices
GetObs = function(df,Node){
  return(df[Node$List,])
}


