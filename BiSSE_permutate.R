library(ape)
library(TreeSim)
library(diversitree)
library(phangorn)
library(hisse)
library(geiger)
library(minqa)
library(R.utils)

gibbs_tree_unpruned <- read.nexus("/Users/fernandovillanea/Documents/Bee_popgen/BiSSE/gibbs_nu.nex")
gibbs_data <- read.csv("/Users/fernandovillanea/Documents/Bee_popgen/BiSSE/gibbs_data.csv", row.names = 1)

name.check(gibbs_tree_unpruned, gibbs_data)

force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}

gibbs_tree<-force.ultrametric(gibbs_tree_unpruned) ## default method
is.ultrametric(gibbs_tree)

#social bees
dropped_bees = c((which(gibbs_data$behav_binary == "unknown")), (which(gibbs_data$behav_binary == "polymorphic")))
kept_bees = c((which(gibbs_data$behav_binary == "eusocial")), (which(gibbs_data$behav_binary == "solitary")))

names = rownames(gibbs_data)
dropped_names = names[dropped_bees]
kept_names = names[kept_bees]

kept_data = gibbs_data[kept_bees,]
names(kept_data) <- kept_names

gibbs_tree_trimmed = drop.tip(gibbs_tree, dropped_names)
plot(gibbs_tree_trimmed)
length(gibbs_tree_trimmed$tip.label)

behavior <-kept_data[,"behav_binary"]
names(behavior) <- rownames(gibbs_data)

is_social <- as.numeric(kept_data == "eusocial")
names(is_social) <- names(kept_data)
is_social

#permutate random pick 15% of states
pval = c()
miss_id_rate = 15/100
for (j in 1:100){
  #reset is_social
  is_social <- as.numeric(kept_data == "eusocial")
  names(is_social) <- names(kept_data)
  
  #choose some tips
  pick = sample(1:(length(gibbs_tree_trimmed$tip.label)), round((length(gibbs_tree_trimmed$tip.label)*miss_id_rate)), replace=FALSE)
  #flip states
  for (i in pick){
    if (is_social[i] == 1) {
      is_social[i] = 0
    }
    else{
      is_social[i] = 1
    }
  }
  
  #run BiSSE and record p-value
  sampling.f = c((length(which(is_social == "0"))/nrow(gibbs_data)),length(which(is_social == "1"))/nrow(gibbs_data))
  nbModel <- make.bisse(gibbs_tree_trimmed, is_social, sampling.f = sampling.f)
  p <- starting.point.bisse(gibbs_tree_trimmed)
  nbModelMLFit <- find.mle(nbModel, p)
  estimated = round(coef(nbModelMLFit))
  cnbModel <- constrain(nbModel, lambda1 ~ lambda0)
  cnbModel <- constrain(cnbModel, mu1 ~ mu0)
  cnbModelMLFit <- find.mle(cnbModel, p[c(-1, -3)], method="subplex")
  anova = anova(nbModelMLFit, constrained = cnbModelMLFit)
  pval = c(pval,anova$"Pr(>|Chi|)"[2])
}
pval_50 = c(pval_50, pval)