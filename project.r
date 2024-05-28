rm(list = ls())

library(igraph)
library(ergm)
library(intergraph)
library(sbm)
library(igraph)
library(ggplot2)


load("workspace.Rdata")
ls()

edgeList = read.table("edge_list.txt", sep = "", head = FALSE)
head(edgeList)
g <- graph_from_edgelist(as.matrix(edgeList), directed=TRUE)

attr = read.table("nodes_attr.txt", sep = "", head = TRUE)
head(attr)


V(g)$continent = attr$Continents
V(g)$gdp = attr$GDP
V(g)$world_partition = attr$World_Partitions
V(g)$xCoordinate = attr$xCoordinates
V(g)$yCoordinate = attr$yCoordinates
V(g)$name = attr$Names

x11(width=100, height = 50)
# Modifica manualmente le posizioni dei nodi
l = matrix(, nrow=length(V(g)), ncol=2)
for (i in 1:length(V(g))) {
  l[i, 1] <- V(g)[i]$xCoordinate* 4 - 2
  l[i, 2] <- 2.25 - V(g)[i]$yCoordinate * 4
}
plot(g, layout = l, rescale=FALSE, vertex.shape = "rectangle", vertex.size=25, vertex.size2= 7)


Y <- as_adjacency_matrix(g)

diag(Y) = NA

rho <- edge_density(g)

reciprocity <- reciprocity(g)

transitivity <- transitivity(g)

odd_rho <- rho / (1 - rho)
odd_transitivity <- transitivity / (1 - transitivity)

tau <- odd_transitivity / odd_rho

assortativity_continent <- assortativity(g, V(g)$continent, directed = TRUE)
assortativity_gdp <- assortativity(g, V(g)$gdp, directed = TRUE)
assortativity_worldpartition <- assortativity(g, V(g)$world_partition, directed = TRUE)

#Nodal Stats

in_degree <- degree(g, mode = "in")
out_degree <- degree(g, mode = "out")

st_in_degree <- degree(g, mode = "in", normalized = TRUE)
st_out_degree <- degree(g, mode = "out", normalized = TRUE)

closeness_in <- closeness(g, mode = "in") #non torna???
closeness_out <- closeness(g, mode = "out")

st_closeness_in <- closeness(g, mode = "in", normalized = TRUE)
st_closeness_out <- closeness(g, mode = "out", normalized = TRUE)

betweenness <- betweenness(g, directed = TRUE)
st_betweenness <- betweenness(g, directed = TRUE, normalized = TRUE)

eigen_centrality <- eigen_centrality(g)$vector
st_eigen_centrality <- eigen_centrality(g, scale = TRUE)$vector

in_centr_degree <- centr_degree(g, loops = FALSE, mode = "in")
out_centr_degree <- centr_degree(g, loops = FALSE, mode = "out")

in_centr_clo <- centr_clo(g, mode = "in")
out_centr_clo <- centr_clo(g, mode = "out")

centr_betw <- centr_betw(g, directed = TRUE)

palette.colors()

summary(in_degree)
hist(in_degree, breaks = length(unique(in_degree)))
in_degree[order(in_degree, decreasing = TRUE)]
summary(st_in_degree)
hist(st_in_degree, breaks = length(unique(st_in_degree)))
plot(g, layout = l, rescale = FALSE, vertex.label.cex = st_in_degree * 5, vertex.size = 0, vertex.label.color = st_in_degree * 8, edge.color = "lightblue")

summary(out_degree)
hist(out_degree, breaks = length(unique(out_degree)))
out_degree[order(out_degree, decreasing = TRUE)]
summary(st_out_degree)
hist(st_out_degree, breaks = length(unique(st_out_degree)))
plot(g, layout = l, rescale = FALSE, vertex.label.cex = st_out_degree * 5, vertex.size = 0, vertex.label.color = st_out_degree * 8, edge.color = "lightblue")

summary(closeness_in)
hist(closeness_in, breaks = length(unique(closeness_in)))
closeness_in[order(closeness_in, decreasing = TRUE)]
summary(st_closeness_in)
hist(st_closeness_in, breaks = length(unique(st_closeness_in)))
plot(g, layout = l, rescale = FALSE, vertex.label.cex = st_closeness_in * 5, vertex.size = 0, vertex.label.color = st_closeness_in * 8, edge.color = "lightblue")

summary(closeness_out)
hist(closeness_out, breaks = length(unique(closeness_out)))
closeness_out[order(closeness_out, decreasing = TRUE)]
summary(st_closeness_out)
hist(st_closeness_out, breaks = length(unique(st_closeness_out)))
plot(g, layout = l, rescale = FALSE, vertex.label.cex = st_closeness_out * 5, vertex.size = 0, vertex.label.color = st_closeness_out * 8, edge.color = "lightblue")

summary(betweenness)
hist(betweenness, breaks = length(unique(betweenness)))
betweenness[order(betweenness, decreasing = TRUE)]
summary(st_betweenness)
hist(st_betweenness, breaks = length(unique(st_betweenness)))
plot(g, layout = l, rescale = FALSE, vertex.label.cex = st_betweenness * 20, vertex.size = 0, vertex.label.color = st_betweenness * 8, edge.color = "lightblue")

summary(eigen_centrality)
hist(eigen_centrality, breaks = length(unique(eigen_centrality)))
eigen_centrality[order(eigen_centrality, decreasing = TRUE)]
summary(st_eigen_centrality)
hist(st_eigen_centrality, breaks = length(unique(st_eigen_centrality)))
plot(g, layout = l, rescale = FALSE, vertex.label.cex = st_eigen_centrality * 3, vertex.size = 0, vertex.label.color = st_eigen_centrality * 8, edge.color = "lightblue")



# SRG

n <- vcount(g)
B <- 1000
m <-  sum(Y, na.rm = TRUE)
#CId.sim = C.sim = c() statistiche da simulare
for (b in 1: B) {
  Y_sim <- matrix(NA, n, n)
  ones <- rep(1, m)
  zeros <- rep(0, n * (n - 1) - m)
  all <- c(ones, zeros)
  Y_sim[col(Y_sim) != row(Y_sim)] <- sample(all, n * (n - 1))
  g_sim <- graph_from_adjacency_matrix(Y_sim)
  #C.sim[b] = transitivity(g.sim)

}

hist(CId.sim, col = "lightgray", main = "Null distribution", xlim = c(0, 0.5))
abline(v = CId.obs, col = "red", lwd=2)
mean(CId.sim >= CId.obs)

### Non-Homogeneous SRG model

aRect_fnc = function(Y, k){

  # Y = adjacency matrix
  # k = n. of steps in the alternating rectangles algorithm
  Y1 <- matrix(c(0,1,1,0), 2, 2)
  Y2 = 1 - Y1

  n = nrow(Y)

  for(s in 1:k){
    # draw 4 distinct indexes
    # two rows and two columns
    ij = sample(1:n,4,replace = F)

    # select the corresponst_in_degreeg sub-matrix
    rows = ij[1:2]
    cols = ij[3:4]
    Yij = Y[rows, cols]

    # perturbation
    if(all(Yij == Y1)) Yij = Y2 else if(all(Yij == Y2))  Yij = Y1

    Y[rows, cols] = Yij
  }

  return(Y)
}

in_sim = out_sim = recip_sim = c()
for(b in 1:1000){
  Y_sim = aRect_fnc(Y, 100)
  g_sim = graph_from_adjacency_matrix(as.matrix(Y_sim))
  cat(sum(Y != Y_sim, na.rm = T), "*", sep="")
  in_sim[b] = sd(degree(g_sim, mode = "in"))
  out_sim[b] = sd(degree(g_sim, mode = "out"))
  recip_sim[b] = reciprocity(g_sim)

}

par(mfrow = c(1,3))
hist(in_sim, breaks = length(unique(in_sim)))
abline(v = sd(in_degree), col = "red", lty = 2)
st_in_sim

hist(out_sim, breaks = length(unique(out_sim)))
abline(v = sd(out_degree), col = "red", lty = 2)
st_out_sim

hist(recip_sim, breaks = length(unique(recip_sim)))
abline(v = reciprocity, col = "red", lty = 2)

# p-value
mean(in_sim > sd(in_degree))
mean(out_sim > sd(out_degree))
mean(recip_sim > reciprocity)




net = network(Y, directed = T)

# net %v% "continent" = V(g)$continent
# net %v% "gdp" = V(g)$gdp
# net %v% "partition" = V(g)$partition
# net %v% "xCoordinate" = V(g)$xCoordinate
# net %v% "yCoordinate" = V(g)$yCoordinate

net %v% "continent" = attr$Continents
net %v% "gdp" = attr$GDP
net %v% "partition" = attr$World_Partitions
net %v% "xCoordinate" = attr$xCoordinates
net %v% "yCoordinate" = attr$yCoordinates
net %v% "name" = attr$Names

mod0 = ergm(net ~ edges, 
            control = control.ergm(seed = 1, checkpoint="mod0/step_%03d.RData"))
summary(mod0)


mod1_sender = ergm(net ~ edges + sender, 
            control = control.ergm(seed = 1, checkpoint="mod1_sender/step_%03d.RData"))
summary(mod1_sender)

mod1_receiver = ergm(net ~ edges + receiver, 
            control = control.ergm(seed = 1, checkpoint="mod1_receiver/step_%03d.RData"))
summary(mod1_receiver)

mod1_mutual = ergm(net ~ edges + mutual, 
              control = control.ergm(seed = 1, checkpoint="mod1_mutual/step_%03d.RData"))
summary(mod1_mutual)

mod2 = ergm(net ~ edges + receiver + mutual,
            verbose = 4,
            control = control.ergm(seed = 1, checkpoint="mod2/step_%03d.RData"))
summary(mod2)
# Non riesco a finalizzare mod2, quindi inutilizzabile per un problema di R

BIC(mod0, mod1_receiver, mod1_mutual)
AIC(mod0, mod1_receiver, mod1_mutual)

# mod1_mutual parrebbe essere il modello migliore per ora

mod3 = ergm(net ~ edges + mutual + 
              nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
              nodefactor("continent") + nodefactor("partition") +
              absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
              nodematch("continent") + nodematch("partition"),
            control = control.ergm(seed = 1, checkpoint="attr_model/step_%03d.RData")) 
summary(mod3)
mcmc.diagnostics(mod3)

BIC(mod1_mutual, mod3)
AIC(mod1_mutual, mod3)

# mod3 sembra il migliore dei due

#Markov
mod4 = ergm(net ~ edges + mutual + 
              nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
              nodefactor("continent") + nodefactor("partition") +
              absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
              nodematch("continent") + nodematch("partition") +
              + istar(2) + ostar(2) + triangle,
            verbose = 4,
            control = control.ergm(seed = 1, checkpoint="markov01/step_%03d.RData")) 
summary(mod4)
mcmc.diagnostics(mod4)
# impossibile da stimare

mod4_noTriangles = ergm(net ~ edges + mutual + 
              nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
              nodefactor("continent") + nodefactor("partition") +
              absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
              nodematch("continent") + nodematch("partition") +
              + istar(2) + ostar(2),
            verbose = 4,
            control = control.ergm(seed = 1, checkpoint="markov02/step_%03d.RData")) 
summary(mod4_noTriangles)
mcmc.diagnostics(mod4_noTriangles)
# impossibile da stimare

mod5 = ergm(net ~ edges + mutual + 
              nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
              nodefactor("continent") + nodefactor("partition") +
              absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
              nodematch("continent") + nodematch("partition") +
              triangle + gwidegree(decay = 1, fixed = TRUE),
            verbose = 4,
            control = control.ergm(seed = 1, checkpoint="markov03/step_%03d.RData")) 
summary(mod5)
mcmc.diagnostics(mod5)
# impossibile da stimare

mod5 = ergm(net ~ edges + mutual + 
              nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
              nodefactor("continent") + nodefactor("partition") +
              absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
              nodematch("continent") + nodematch("partition") +
              triangle + gwodegree(decay = 1, fixed = TRUE),
            verbose = 4,
            control = control.ergm(seed = 1, checkpoint="markov03gw01/step_%03d.RData")) 
summary(mod5)
mcmc.diagnostics(mod5)
# impossibile da stimare

mod5 = ergm(net ~ edges + mutual + 
              nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
              nodefactor("continent") + nodefactor("partition") +
              absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
              nodematch("continent") + nodematch("partition") +
              triangle + gwidegree(decay = 1, fixed = TRUE) + gwodegree(decay = 1, fixed = TRUE),
            verbose = 4,
            control = control.ergm(seed = 1, checkpoint="markov03gw02/step_%03d.RData")) 
summary(mod5)
mcmc.diagnostics(mod5)
# impossibile da stimare


mod6 = ergm(net ~ edges + mutual + 
              nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
              nodefactor("continent") + nodefactor("partition") +
              absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
              nodematch("continent") + nodematch("partition") +
              gwesp(decay = 1, fixed = T) + gwdsp(decay = 1, fixed = T),
            verbose = 4,
            control = control.ergm(seed = 1, checkpoint="mod6/step_%03d.RData")) 
summary(mod6)
mcmc.diagnostics(mod6)
# impossibile da stimare

mod6_gwdsp = ergm(net ~ edges + mutual + 
              nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
              nodefactor("continent") + nodefactor("partition") +
              absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
              nodematch("continent") + nodematch("partition") +
              gwdsp(decay = 1, fixed = T),
            verbose = 4,
            control = control.ergm(seed = 1, checkpoint="mod6_u/step_%03d.RData")) 
summary(mod6_gwdsp)
mcmc.diagnostics(mod6_gwdsp)

# Devo provare mod6 solo con gwesp pe capire

mod6_gwesp = ergm(net ~ edges + mutual + 
                nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
                nodefactor("continent") + nodefactor("partition") +
                absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
                nodematch("continent") + nodematch("partition") +
                gwesp(decay = 1, fixed = T),
              verbose = 4,
              control = control.ergm(seed = 1, checkpoint="mod6_gwesp/step_%03d.RData")) 
summary(mod6_gwesp)
mcmc.diagnostics(mod6_gwesp)


sim = simulate(mod6_u, nsim = 100, verbose = TRUE, seed = 1)

fnc = function(xx){
  ig = asIgraph(xx)
  tr = transitivity(ig)
  ideg = sd(degree(ig, mode = "in"))
  odeg = sd(degree(ig, mode = "out"))
  return(c(tr, ideg, odeg))
}

null.distr = matrix(,100,3)
for(b in 1:100){
  null.distr[b,]  = fnc(sim[[b]])
}
dev.new()
par(mfrow = c(3,1))
hist(unlist(null.distr[,1]), xlab = "transitivity"); abline(v = transitivity(g), col = "red")
hist(unlist(null.distr[,2]), xlab = "in-degree"); abline(v = sd(degree(g, mode = "in")), col = "red")
hist(unlist(null.distr[,3]), xlab = "out-degree"); abline(v = sd(degree(g, mode = "out")), col = "red")

mean(transitivity(g) > null.distr[, 1])
mean(sd(degree(g, mode = "in")) > null.distr[, 2])
mean(sd(degree(g, mode = "out")) > null.distr[, 3])

# SBM
plotMyMatrix(as.matrix(Y), dimLabels = list(row = 'nations', col = 'nations'))

sbm1 = estimateSimpleSBM(as.matrix(Y), "bernoulli", dimLabels = 'nations', 
                         estimOptions = list(verbosity = 1))

# let us look at the results
sbm1

# selected number of blocks
sbm1$nbBlocks
# prior block probabilities
sbm1$blockProp
# connectivity parameters
round(sbm1$connectParam$mean,3)

# a clearer representation
# ------------------------
plot(sbm1, type = "data")
# nodes are ordered wrt to the block they belong to and blocks are highlighted
plot(sbm1, type = "expected")
# fitted connection probabilities
plot(sbm1, type = "meso")
# fitted connection probabilities 

# info on all estimated model is given in 
sbm1$storedModels

sbm1$memberships

plot(g, vertex.color = sbm1$memberships)
