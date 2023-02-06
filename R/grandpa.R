#' @import tidyr
#' @import dplyr 
#' @import igraph
#' @import ggplot2
NULL

#' Conducts a GRANDPA procedure
#' 
#' Uses the GRANDPA framework to generate a network with matching attribute probabilities
#' @param G an igraph object containing the original graph. Attributes which will be in the modelling procedure must contain the word 'Label'
#' @param nt the number of desired nodes in the generated graph
#' @param mt the number of desired edges in the generated graph
#' @param preventSelf logical; should self connections be prevented?
#' @param preventDups logical; should duplicate connections be prevented?
#' @param augment logical; should degree augmentation occur?
#' @param augmentCommunity logical; should communited be augmented?
#' @param nbins the number of degree bins used in the degree augmentation
#' @param degType string; the type of degree to be used in the degree augmentation. See ?degree for details
#' @param seed a random seed set for reproducibility
#' 
#' @return an igraph object with attribute and augmentation labels
#' @examples 
#' 
#' G = make_graph("Zachary")
#' V(G)$degree<-degree(G,mode="all")
#' 
#' V(G)$Label1<-rep("A",vcount(G))
#' V(G)$Label1[V(G)$degree>median(V(G)$degree)]<-"B"
#' 
#' GSim<-grandpa(G)
#' 
#' @export
grandpa<-function(G,nt=vcount(G),mt=ecount(G),preventSelf=T,preventDups=T,augment=F,augmentCommunity=F,nbins=3,degType="all",seed=0){
  
  if(is_igraph(G)==F){
    stop("Graph is not an igraph object")
  }
  
  #Add augmentation label
  if(augment){
    G<-addAugment(G,nbins,degType)
  }
  
  if(augmentCommunity){
    G<-addCommunityAugment(G)
  }
  
  #check if a variable named "label" exists here
  if(any(grepl("label",vertex_attr_names(CMSSub),ignore.case=T))==F){
    stop("At least one vertex attribute must contain the word label")
  }
  
  #calculate observed probabilities of attributes for vertices and edges
  obsProb<-calcObsProbs(G)
  observedLabels<-obsProb[[1]]
  observedEdgeLabels<-obsProb[[2]]
  
  #sample vertices
  sampled_v<-sampleVertices(observedLabels,nt,seed)
  
  #initialize graph
  Gt<-initializeTargetG(sampled_v,nt,d=is_directed(G))
  
  #map category labels to vertices
  C2V<-createMap(sampled_v)
  
  #sample edges
  sampled_e<-sampleEdges(observedEdgeLabels,mt,seed)
  
  #put vertices and edges into graph
  Gt<-simGraph(Gt,sampled_e,C2V,seed,preventSelf,preventDups)
  
  return(Gt)
}

#' Calculates the observed edge and vertex probabilities
#' 
#' This function is typically called by the grandpa function rather than directly.
#' @param G An igraph object, containing labels indicated with the word 'Label' that should be used to calculate edge and vertex probabilities
#' 
#' @return a list of tidy tables containing the vertex probabilities and edge probabilities respectively
#' @examples 
#' 
#' G = make_graph("Zachary")
#' V(G)$degree<-degree(G,mode="all")
#' 
#' V(G)$Label1<-rep("A",vcount(G))
#' V(G)$Label1[V(G)$degree>median(V(G)$degree)]<-"B"
#' 
#' obsProbTabs<-calcObsProbs(G)
#' vertexProbs<-obsProbTabs[[1]]
#' edgeProbs<-obsProbTabs[[2]]
#' @export
calcObsProbs<-function(G){
  #This function calculates the observed probabilities of combinations of attributes occurring on vertices, as well as edges connecting those vertices
  
  #get all Label variables
  Labels<-list.vertex.attributes(G)[grepl("label",list.vertex.attributes(G),ignore.case = T)]
  
  #find all combinations
  observedLabels<-vertex_attr(G, Labels[1])
  
  if(length(Labels)>1){
    for(i in 2:length(Labels)){
      observedLabels<-cbind(observedLabels,vertex_attr(G, Labels[i]))
    }}
  
  observedLabels<-data.frame(observedLabels)
  colnames(observedLabels)<-Labels
  rownames(observedLabels)<-V(G)$name
  #count combos and find probabilities
  observedLabels_summary = observedLabels %>% count(across()) %>% mutate(p=n/sum(n)) %>% arrange(desc(p))
  
  #get all edges
  Es<-as_edgelist(G,names=F)
  
  #make combos, count, find probability
  if(length(Labels)==1){
    observedEdgeLabels<-data.frame(observedLabels[Es[,1],],observedLabels[Es[,2],]) %>%
      count(across()) %>% mutate(p=n/sum(n)) %>% arrange(desc(p))
  } else {
    observedEdgeLabels<-data.frame(do.call(paste,c(observedLabels,sep=""))[Es[,1]],
                                   do.call(paste,c(observedLabels,sep=""))[Es[,2]]) %>% 
      count(across()) %>% mutate(p=n/sum(n)) %>% arrange(desc(p))
  }
  
  
  return(list(observedLabels_summary,observedEdgeLabels))
}

#' Creates the vertices that will be used in the generated graph
#' 
#' This function is typically called by the grandpa function rather than directly.
#' @param observedLabels a dataframe containing the observed  label combinations and their probabilities from the original graph
#' @param nt the number of desired nodes in the generated graph
#' @param seed the random seed to be used for reproducibility
#' 
#' @return a dataframe containing node attributes
#' @examples 
#' 
#' G = make_graph("Zachary")
#' V(G)$degree<-degree(G,mode="all")
#' 
#' V(G)$Label1<-rep("A",vcount(G))
#' V(G)$Label1[V(G)$degree>median(V(G)$degree)]<-"B"
#' 
#' obsProbTabs<-calcObsProbs(G)
#' vertexProbs<-obsProbTabs[[1]]
#' edgeProbs<-obsProbTabs[[2]]
#' 
#' sampledV<-sampleVertices(vertexProbs,50,0)
#' @export
sampleVertices<-function(observedLabels,nt,seed){
  set.seed(seed)
  observedLabels$addToBag<-round(observedLabels$p*nt,0)
  sampled_v<-observedLabels %>% select((contains("label",ignore.case=T)),addToBag) %>% uncount(addToBag)
  sampled_v<-data.frame(sampled_v[sample(1:nt),])
  colnames(sampled_v)<-colnames(observedLabels)[grepl("label",colnames(observedLabels),ignore.case = T)]
  sampled_v$ix<-1:nt
  return(sampled_v)
}

#' Intializes an empty graph with the generated vertices
#' 
#' This function is typically called by the grandpa function rather than directly.
#' @param sampled_v a dataframe containing the sampled vertices and their corresponding labels on each row
#' @param nt the number of desired nodes in the generated graph
#' @param d logical; is the graph directed?
#' 
#' @return an empty igraph object with the desired vertices
#' @examples 
#' 
#' G = make_graph("Zachary")
#' V(G)$degree<-degree(G,mode="all")
#' 
#' V(G)$Label1<-rep("A",vcount(G))
#' V(G)$Label1[V(G)$degree>median(V(G)$degree)]<-"B"
#' 
#' obsProbTabs<-calcObsProbs(G)
#' vertexProbs<-obsProbTabs[[1]]
#' edgeProbs<-obsProbTabs[[2]]
#' 
#' sampledV<-sampleVertices(vertexProbs,50,0)
#' initG<-initializeTargetG(sampledV,50,FALSE)
#' @export
initializeTargetG<-function(sampled_v,nt,d){
  Gt<-make_empty_graph(nt,directed=d)
  
  #set attributes
  for(cn in colnames(sampled_v)[-ncol(sampled_v)]){
    Gt=set_vertex_attr(Gt,cn,1:nt,value=as.character(sampled_v[,cn]))
  }
  
  if(length(colnames(sampled_v))>2){
    V(Gt)$SummaryLabel<-do.call(paste,c(sampled_v[,1:ncol(sampled_v)-1],sep=""))
  } else {Gt=set_vertex_attr(Gt,"SummaryLabel",1:nt,value=as.character(sampled_v[,cn]))}
  
  return(Gt)
}

#' Creates the map which lists vertix IDs based on attribute labels
#' 
#' This function is typically called by the grandpa function rather than directly.
#' @param sampled_v a dataframe containing the sampled vertices and their corresponding labels on each row
#' 
#' @return a list which contains a mapping of all vertex IDs based on label combinations
#' @examples 
#' 
#' G = make_graph("Zachary")
#' V(G)$degree<-degree(G,mode="all")
#' 
#' V(G)$Label1<-rep("A",vcount(G))
#' V(G)$Label1[V(G)$degree>median(V(G)$degree)]<-"B"
#' 
#' obsProbTabs<-calcObsProbs(G)
#' vertexProbs<-obsProbTabs[[1]]
#' edgeProbs<-obsProbTabs[[2]]
#' 
#' sampledV<-sampleVertices(vertexProbs,50,0)
#' map<-createMap(sampledV)
#' @export
createMap<-function(sampled_v){
  #Create mapping
  C2V<-suppressMessages(sampled_v %>% group_by(across(!c(ix))) %>% summarise(map=list(ix)))
  names(C2V$map)<-do.call(paste,c(C2V[,-ncol(C2V)],sep=""))
  return(C2V)
}

#' Creates the edges that will be used in the generated graph
#' 
#' This function is typically called by the grandpa function rather than directly.
#' @param observedEdgeLabels a dataframe containing the observed edge label combinations and their probabilities from the original graph
#' @param mt the number of desired edges in the generated graph
#' @param seed the random seed to be used for reproducibility
#' 
#' @return a dataframe containing edges and their nodel attributes on either end
#' @examples 
#' 
#' G = make_graph("Zachary")
#' V(G)$degree<-degree(G,mode="all")
#' 
#' V(G)$Label1<-rep("A",vcount(G))
#' V(G)$Label1[V(G)$degree>median(V(G)$degree)]<-"B"
#' 
#' obsProbTabs<-calcObsProbs(G)
#' vertexProbs<-obsProbTabs[[1]]
#' edgeProbs<-obsProbTabs[[2]]
#' 
#' sampledE<-sampleVertices(edgeProbs,75,0)
#' @export
sampleEdges<-function(observedEdgeLabels,mt,seed){
  #draw edges
  observedEdgeLabels$addToBag<-round(observedEdgeLabels$p*mt,0)
  sampled_e<-observedEdgeLabels %>% select(1:2,addToBag) %>% uncount(addToBag)
  sampled_e<-sampled_e[sample(1:mt),]
  sampled_e$ix<-1:mt
  colnames(sampled_e)<-paste("Label",1:2,sep="")
  return(sampled_e)
}

#' Procedure to add edges to the generated graph
#' 
#' This function is typically called by the grandpa function rather than directly.
#' @param Gt The empty graph
#' @param sampled_e a dataframe containing the edges for the generated graph and their node attributes
#' @param C2V a mapping of the nodal attribute combinations and vertex IDs
#' @param seed the random seed to be used for reproducibility
#' @param preventSelf logical; should self connections be prevented?
#' @param preventDups logical; should duplicate connections be prevented
#' 
#' @return a generated igraph object
#' @examples 
#' 
#' G = make_graph("Zachary")
#' V(G)$degree<-degree(G,mode="all")
#' 
#' V(G)$Label1<-rep("A",vcount(G))
#' V(G)$Label1[V(G)$degree>median(V(G)$degree)]<-"B"
#' 
#' obsProbTabs<-calcObsProbs(G)
#' vertexProbs<-obsProbTabs[[1]]
#' edgeProbs<-obsProbTabs[[2]]
#' 
#' sampledV<-sampleVertices(vertexProbs,50,0)
#' initG<-initializeTargetG(sampledV,50,FALSE)
#' map<-createMap(sampledV)
#' sampledE<-sampleVertices(edgeProbs,75,0)
#' 
#' GSim<-simGraph(initG,sampledE,map,0,T,T)
#' @export
simGraph<-function(Gt,sampled_e,C2V,seed,preventSelf,preventDups){
  set.seed(seed)
  first=T
  
  for(i in 1:nrow(sampled_e)){
    lab1<-as.character(sampled_e[i,1])
    lab2<-as.character(sampled_e[i,2])
    
    if(length(C2V$map[[lab1]])>1){
      v1<-sample(C2V$map[[lab1]],1,prob=1/(degree(Gt,v=C2V$map[[lab1]])+1e-4))
    } else {v1=C2V$map[[lab1]]}
    
    if(length(C2V$map[[lab2]])>1){    #probability inversely proportional to degree in target graph
      v2<-sample(C2V$map[[lab2]],1,prob=1/(degree(Gt,v=C2V$map[[lab2]])+1e-4))
    } else {v2=C2V$map[[lab2]]}
    
    #prevent self connecting for now, add as option to override later
    if(preventSelf&v1==v2){
      v_ix<-preventSelfConnections(v1,v2,lab1,lab2,C2V)
      v1<-v_ix[1]
      v2<-v_ix[2]
    }
    
    #avoid drawing the same edge twice
    if(length(C2V$map[[lab1]])==1 & length(C2V$map[[lab2]]==1) & are.connected(Gt, v1, v2) & preventDups){
      if(first==T){
        warning("Could not prevent duplicate connections. Consider simplifying to unique edges")
        first=F
      }
      next
    }
    
    if(preventDups&are.connected(Gt, v1, v2)){
      v_ix<-preventDuplicateConnections(v1,v2,lab1,lab2,C2V,Gt)
      v1<-v_ix[1]
      v2<-v_ix[2]
    }
    
    #add edge in igraph
    Gt<-add_edges(Gt,c(v1,v2))
  }
  return(Gt)
}

#' Add degree augmentation label
#' 
#' This function is typically called by the grandpa function rather than directly.
#' @param G an igraph object
#' @param nbins the number of degree bins used in the degree augmentation
#' @param degType string; the type of degree to be used in the degree augmentation. See ?degree for details
#' 
#' @return an igraph object containing an augmentation label corresponding to degree
#' @examples 
#' 
#' G <- make_graph("Zachary")
#' G <- addAugment(G,3,"all")
#' @export
addAugment<-function(G,nbins,degType){
  #currently defaulting to all, can append later
  deg<-degree(G,mode=degType)
  V(G)$augLabel<-ggplot2::cut_interval(deg,nbins)
  return(G)
}

#' Add community augmentation label
#' 
#' This function is typically called by the grandpa function rather than directly. By default, this function implements louvain clustering to detect community membership. Note, it is recommended to fit a community algorithm outside of the GRANDPA framework and add a vertext label called 'CommunityLabel' for best results.
#' @param G an igraph object
#' 
#' @return an igraph object containing an augmentation label corresponding to community membership
#' @examples 
#' 
#' G <- make_graph("Zachary")
#' G <- addCommunityAugment(G)
#' @export
addCommunityAugment<-function(G){
  if(is_directed(G)){
    warning("Communities are detected on an undirected graph")
    tempG<-as.undirected(G,mode="collapse")
  } else {tempG<-G}
  Gclust<-cluster_louvain(tempG,resolution=resolution)
  V(G)$CommunityLabel<-membership(Gclust)
  return(G)
}

#' Prevents self connections during graph generation
#' 
#' This function is typically called by the grandpa function and should not be used directly.
#' @param v1 a vertex id
#' @param v2 a second vertex id
#' @param lab1 the label corresponding to vertex id 1
#' @param lab2 the label corresponding to vertex id 2
#' @param C2V a mapping of label categories to vertex IDs
#' 
#' @return updated vertex IDs
#' @export
preventSelfConnections<-function(v1,v2,lab1,lab2,C2V){
  while(v1==v2&length(C2V$map[[lab1]])>1){
    v1<-sample(C2V$map[[lab1]],1)
  }
  while(v1==v2&length(C2V$map[[lab2]])>1){
    v2<-sample(C2V$map[[lab2]],1)
  }
  if(v1==v2){
    warning("Could not prevent self connection")
  }
  return(c(v1,v2))
}

#' Prevents duplicate connections during graph generation
#' 
#' This function is typically called by the grandpa function and should not be used directly.
#' @param v1 a vertex id
#' @param v2 a second vertex id
#' @param lab1 the label corresponding to vertex id 1
#' @param lab2 the label corresponding to vertex id 2
#' @param C2V a mapping of label categories to vertex IDs
#' @param Gt a graph being generated
#' 
#' @return updated vertex IDs
#' @export
preventDuplicateConnections<-function(v1,v2,lab1,lab2,C2V,Gt){
  df<-expand_grid(C2V$map[[lab1]],C2V$map[[lab2]])
  colnames(df)<-c("V1","V2")
  df<-df %>% rowwise() %>% mutate(areConnected=are.connected(Gt, V1, V2)) %>% filter(areConnected==F)
  
  if(nrow(df)==F){
    warning("Could not prevent duplicate connections. Consider simplifying to unique edges")
    return(c(v1,v2))
  }
  else{
    ix<-sample(1:nrow(df),1)
    return(as.numeric(df[ix,1:2]))
  }
}

#' Calculates the error between two degree distributions of graphs
#' 
#' This function takes an original graph and a simulated graph and calculates the normalized root mean square error between their degree distributions.
#' @param G the original graph
#' @param Gt the generated graph
#' 
#' @return numerical; the NRMSE between the degree distributions
#' @export
degree_error<-function(Gs,Gt){
  
  sourceCCDF<-degree_distribution(G,cumulative=T)
  targetCCDF<-degree_distribution(Gt,cumulative=T)
  
  source_f<-approxfun(sourceCCDF)
  target_f<-approxfun(targetCCDF)
  
  overlap_x<-c(1,min(length(sourceCCDF),length(targetCCDF)))
  i<-seq(1,max(overlap_x),length.out=500)
  
  h1<-source_f(i)-target_f(i) 
  error<-sqrt(sum(h1^2)/length(h1))/mean(target_f(i))
  return(error)
}