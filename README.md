# GeneRAting Networks using Degree and Property Augmentation

GRANDPA is a graph generating algorithm which uses observed attribute frequencies and topological features to simulate graphs from real world network data.

To install, use the following:

```{r}
devtools::install_github("https://github.com/CarlyBobak/grandpa")
library(grandpa)
library(igraph)
```

Generate a base graph:

```{r}
G = make_graph("Zachary")
coords = layout_with_fr(G)
# plot the graph
plot(G, layout=coords, vertex.label=NA, vertex.size=10)
```

Add associated vertex attributes using the following:

```{r}
#Label possibilities
L1<-c(1,2,3)
L2<-c("A","B","C","D")
L<-list(L1,L2)
L

#calculate degree and centrality
V(G)$degree<-degree(G,mode="all")
V(G)$cent<-unlist(centr_clo(G,mode="total")$res) #closeness centrality

Qcent<-quantile(V(G)$cent)
Qdeg<-quantile(V(G)$degree)

lab1<-c()
lab2<-c()

for(i in V(G)){
  #For L1
  #Consider: a more flexible way of defining degree cutoffs for different graphs
  if(V(G)$degree[i]<Qdeg[3]){
    lab1[i]<-sample(L[[1]],1,prob=c(0.8,0.1,0.1))
  } else if(V(G)$degree[i]<Qdeg[4]){
    lab1[i]<-sample(L[[1]],1,prob=c(0.4,0.3,0.3))
  } else{
    lab1[i]<-sample(L[[1]],1,prob=c(0.1,0.2,0.7))
  }
  
  #For L2
  if(V(G)$cent[i]<Qcent[2]){
    lab2[i]<-sample(L[[2]],1,prob=c(0.25,0.25,0.25,0.25))
  } else if(V(G)$cent[i]<Qcent[3]){
    lab2[i]<-sample(L[[2]],1,prob=c(0.25,0.25,0.4,0.1))
  } else if(V(G)$cent[i]<Qcent[4]){
    lab2[i]<-sample(L[[2]],1,prob=c(0.25,0.5,0.15,0.1))
  } else{
    lab2[i]<-sample(L[[2]],1,prob=c(0.75,0.1,0.1,0.05))
  }
}

V(G)$Label1<-lab1
V(G)$Label2<-lab2

#plot with label information
V(G)$size<-as.numeric(as.factor(lab2))*3
V(G)$color<-as.factor(lab1)

plot(G, layout=coords,vertex.label.cex=0.5)
```

Obtain clusters:
```{r}
c = cluster_edge_betweenness(G)
V(G)$CommunityLabel<-membership(c)
```

Use GRANDPA to generate a new graph:
```{r}
GSim<-grandpa(G,augment = T,nbins=10)
GSim<-simplify(GSim)
V(GSim)$size<-as.numeric(as.factor(V(GSim)$Label2))*3
V(GSim)$color<-as.factor(V(GSim)$Label1)
plot(GSim, layout = layout_with_fr(GSim), vertex.label.cex=0.5,mark.groups = cluster_edge_betweenness(GSim))
```
