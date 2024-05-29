rm(list = ls())

library(igraph)
library(ergm)
library(intergraph)
library(sbm)
library(ggplot2)
library(gridExtra)

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
continent_colors <- palette.colors(9)[-1]
plot(g, layout = l, rescale=FALSE, vertex.shape = "rectangle", vertex.size=19, vertex.size2= 3, vertex.label.font = 2, vertex.color = continent_colors[V(g)$continent])
continent_names <- c("Africa", "Asia", "Europe", "North America", "Oceania", "South America")
legend("topright", legend = continent_names, fill = continent_colors, cex = 2)

plot(g, layout = l, rescale=FALSE, vertex.size=V(g)$gdp/2000, vertex.label.font = 2)

world_partition_colors <- palette.colors(9)[-1]
plot(g, layout = l, rescale=FALSE, vertex.label.font = 2, vertex.color = world_partition_colors[V(g)$world_partition], vertex.size=10)
world_partition_names <- c("Core", "Semiperiphery", "Periphery")
legend("topright", legend = world_partition_names, fill = world_partition_colors, cex = 2)


Y <- as_adjacency_matrix(g)

diag(Y) <- NA

rho <- edge_density(g)

reciprocity <- reciprocity(g)

transitivity <- transitivity(g)

odd_rho <- rho / (1 - rho)
odd_transitivity <- transitivity / (1 - transitivity)

tau <- odd_transitivity / odd_rho

assortativity_continent <- assortativity(g, V(g)$continent, directed = TRUE)
assortativity_gdp <- assortativity(g, V(g)$gdp, directed = TRUE)
assortativity_worldpartition <- assortativity(g, V(g)$world_partition, directed = TRUE)

# Create a table with names and values
table <- data.frame(
  Names = c("Edge Density", "Reciprocity", "Transitivity", "Odd-Ratio of Edge Density", "Odd-Ratio of Transitivity", "Tau", "Assortativity (Continent)", "Assortativity (GDP)", "Assortativity (World Partition)"),
  Values = c(rho, reciprocity, transitivity, odd_rho, odd_transitivity, tau, assortativity_continent, assortativity_gdp, assortativity_worldpartition)
)
# Save the table as an image
png("table.png", width = 800*3, height = 400*3, res=72*3)
grid.table(table)
dev.off()


#Nodal Stats

in_degree <- degree(g, mode = "in")
out_degree <- degree(g, mode = "out")
st_in_degree <- degree(g, mode = "in", normalized = TRUE)
st_out_degree <- degree(g, mode = "out", normalized = TRUE)

closeness_in <- closeness(g, mode = "in")
closeness_out <- closeness(g, mode = "out")
st_closeness_in <- closeness(g, mode = "in", normalized = TRUE)
st_closeness_out <- closeness(g, mode = "out", normalized = TRUE)

betweenness <- betweenness(g, directed = TRUE)
st_betweenness <- betweenness(g, directed = TRUE, normalized = TRUE)
eigen_centrality <- eigen_centrality(g)$vector

in_centr_degree <- centr_degree(g, loops = FALSE, mode = "in")
out_centr_degree <- centr_degree(g, loops = FALSE, mode = "out")
in_centr_clo <- centr_clo(g, mode = "in")
out_centr_clo <- centr_clo(g, mode = "out")

centr_betw <- centr_betw(g, directed = TRUE)

palette.colors()

summary(in_degree)
in_degree[order(in_degree, decreasing = TRUE)]
summary(st_in_degree)

top_nodes <- names(sort(st_in_degree, decreasing = TRUE)[1:10])
label_colors <- rep("black", length(V(g)))
label_colors[V(g)$name %in% top_nodes] <- "red"
plot(g, layout = l, rescale = FALSE, vertex.label.cex = st_in_degree * 5, vertex.size = 0, vertex.label.color = label_colors, edge.color = "lightblue")

summary(out_degree)
out_degree[order(out_degree, decreasing = TRUE)]
summary(st_out_degree)

top_nodes_out <- names(sort(st_out_degree, decreasing = TRUE)[1:6])
label_colors_out <- rep("black", length(V(g)))
label_colors_out[V(g)$name %in% top_nodes_out] <- "red"
plot(g, layout = l, rescale = FALSE, vertex.label.cex = st_out_degree * 5, vertex.size = 0, vertex.label.color = label_colors_out, edge.color = "lightblue")

summary(closeness_in)
closeness_in[order(closeness_in, decreasing = TRUE)]
summary(st_closeness_in)

top_nodes <- names(sort(st_closeness_in, decreasing = TRUE)[1:9])
label_colors <- rep("black", length(V(g)))
label_colors[V(g)$name %in% top_nodes] <- "red"
plot(g, layout = l, rescale = FALSE, vertex.label.cex = st_closeness_in * 4, vertex.size = 0, vertex.label.color = label_colors, edge.color = "lightblue")

summary(closeness_out)
closeness_out[order(closeness_out, decreasing = TRUE)]
summary(st_closeness_out)

top_nodes_out <- names(sort(st_closeness_out, decreasing = TRUE)[1:5])
label_colors_out <- rep("black", length(V(g)))
label_colors_out[V(g)$name %in% top_nodes_out] <- "red"
plot(g, layout = l, rescale = FALSE, vertex.label.cex = st_closeness_out * 3, vertex.size = 0, vertex.label.color = label_colors_out, edge.color = "lightblue")

summary(betweenness)
betweenness[order(betweenness, decreasing = TRUE)]
summary(st_betweenness)

top_nodes <- names(sort(st_betweenness, decreasing = TRUE)[1:5])
label_colors <- rep("black", length(V(g)))
label_colors[V(g)$name %in% top_nodes] <- "red"
plot(g, layout = l, rescale = FALSE, vertex.label.cex = st_betweenness * 20, vertex.size = 0, vertex.label.color = label_colors, edge.color = "lightblue")

summary(eigen_centrality)
eigen_centrality[order(eigen_centrality, decreasing = TRUE)]

top_nodes <- names(sort(eigen_centrality, decreasing = TRUE)[1:5])
label_colors <- rep("black", length(V(g)))
label_colors[V(g)$name %in% top_nodes] <- "red"
plot(g, layout = l, rescale = FALSE, vertex.label.cex = st_eigen_centrality * 3, vertex.size = 0, vertex.label.color = label_colors, edge.color = "lightblue")


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

#SRG

mod0 = ergm(net ~ edges,
            control = control.ergm(seed = 1, checkpoint="mod0/step_%03d.RData"))
summary(mod0)

#NH-SRG

mod1_receiverAndSender = ergm(net ~ edges + receiver + sender,
                     control = control.ergm(seed = 1))
summary(mod1_receiverAndSender)

mod1_sender = ergm(net ~ edges + sender,
            control = control.ergm(seed = 1, checkpoint="mod1_sender/step_%03d.RData"))
summary(mod1_sender)

mod1_receiver = ergm(net ~ edges + receiver,
            control = control.ergm(seed = 1, checkpoint="mod1_receiver/step_%03d.RData"))
summary(mod1_receiver)

#p1 model

mod2 = ergm(net ~ edges + receiver + mutual,
            verbose = 4,
            control = control.ergm(seed = 1, checkpoint="mod2/step_%03d.RData"))
summary(mod2)
# Non riesco a finalizzare mod2, quindi inutilizzabile per un problema di R (neanche da fisso)

mod2_mutual = ergm(net ~ edges + mutual,
              control = control.ergm(seed = 1))
summary(mod2_mutual)

BIC(mod0, mod1_receiver, mod2_mutual)
AIC(mod0, mod1_receiver, mod2_mutual)

# mod2_mutual parrebbe essere il modello migliore per ora

mod3 = ergm(net ~ edges + mutual +
              nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
              nodefactor("continent") + nodefactor("partition") +
              absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
              nodematch("continent") + nodematch("partition"),
            control = control.ergm(seed = 1, checkpoint="attr_model/step_%03d.RData"))
summary(mod3)
mcmc.diagnostics(mod3)

BIC(mod0, mod2_mutual, mod3)
AIC(mod0, mod2_mutual, mod3)

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
# impossibile da stimare (anche da fisso)

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
# impossibile da stimare (anche da fisso) degenera

mod5 = ergm(net ~ edges + mutual +
              nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
              nodefactor("continent") + nodefactor("partition") +
              absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
              nodematch("continent") + nodematch("partition") +
              gwidegree(decay = 1, fixed = TRUE) + gwodegree(decay = 1, fixed = TRUE),
            verbose = 4,
            control = control.ergm(seed = 1, checkpoint="markov03/step_%03d.RData"))
summary(mod5)
mcmc.diagnostics(mod5)
# impossibile da stimare

mod5_onlyGwod = ergm(net ~ edges + mutual +
              nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
              nodefactor("continent") + nodefactor("partition") +
              absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
              nodematch("continent") + nodematch("partition") +
              gwodegree(decay = 1, fixed = TRUE),
            verbose = 4,
            control = control.ergm(seed = 1, checkpoint="markov03gw01/step_%03d.RData"))
summary(mod5_onlyGwod)
mcmc.diagnostics(mod5_onlyGwod)
# impossibile da stimare

mod5_onlyGwid = ergm(net ~ edges + mutual +
              nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
              nodefactor("continent") + nodefactor("partition") +
              absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
              nodematch("continent") + nodematch("partition") +
              gwidegree(decay = 1, fixed = TRUE),
            verbose = 4,
            control = control.ergm(seed = 1))
summary(mod5_onlyGwid)
mcmc.diagnostics(mod5_onlyGwid)

mod5_onlyGwid_2 = ergm(net ~ edges + mutual +
                       nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
                       nodefactor("continent") + nodefactor("partition") +
                       absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
                       nodematch("continent") +
                       gwidegree(decay = 1, fixed = TRUE),
                     verbose = 4,
                     control = control.ergm(seed = 1, checkpoint="markov03gw03/step_%03d.RData"))
summary(mod5_onlyGwid_2)
mcmc.diagnostics(mod5_onlyGwid_2)

AIC(mod3, mod5_onlyGwid, mod5_onlyGwid_2)
BIC(mod3, mod5_onlyGwid, mod5_onlyGwid_2)

# mod5_onlyGwid_2 sembra migliore

#social circuit model

mod6 = ergm(net ~ edges + mutual +
              nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
              nodefactor("continent") + nodefactor("partition") +
              absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
              nodematch("continent") + nodematch("partition") +
              gwidegree(decay = 1, fixed = TRUE) +
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
            control = control.ergm(seed = 1, checkpoint="mod6_gwdsp/step_%03d.RData"))
summary(mod6_gwdsp)
mcmc.diagnostics(mod6_gwdsp)

mod6_gwesp = ergm(net ~ edges + mutual +
                nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
                nodefactor("continent") + nodefactor("partition") +
                absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
                nodematch("continent") + nodematch("partition") +
                gwesp(decay = 1, fixed = T),
              verbose = 4,
              control = control.ergm(seed = 1, checkpoint="mod6_gwesp/step_%03d.RData", resume = "mod6_gwesp/step_026.RData"))
summary(mod6_gwesp)
mcmc.diagnostics(mod6_gwesp)
# impossibile da stimare

mod6_gwdsp_2 = ergm(net ~ mutual +
                      nodecov("gdp") + nodecov("xCoordinate") + nodecov("yCoordinate") +
                      nodefactor("continent") + nodefactor("partition") +
                      absdiff("gdp") + absdiff("xCoordinate") + absdiff("yCoordinate") +
                      nodematch("continent") + nodematch("partition") +
                      gwdsp(decay = 1, fixed = T),
                    verbose = 4,
                    control = control.ergm(seed = 1, checkpoint="mod6_gwdsp_2/step_%03d.RData"))
summary(mod6_gwdsp_2)
mcmc.diagnostics(mod6_gwdsp_2)

AIC(mod5_onlyGwid_2, mod6_gwdsp, mod6_gwdsp_2)
BIC(mod5_onlyGwid_2, mod6_gwdsp, mod6_gwdsp_2)

# mod6_gwdsp_2 sembra essere il migliore anche se di veramente poco

fnc = function(xx){
  ig = asIgraph(xx)
  tr = transitivity(ig)
  ideg = sd(degree(ig, mode = "in"))
  odeg = sd(degree(ig, mode = "out"))
  dens = edge_density(ig)
  return(c(tr, ideg, odeg, dens))
}

sim2 = simulate(mod6_gwdsp_2, nsim = 100, verbose = TRUE, seed = 1)

null.distr = matrix(,100,4)
for(b in 1:100){
  null.distr[b,]  = fnc(sim2[[b]])
}
dev.new()
par(mfrow = c(4,1))
hist(unlist(null.distr[,1]), xlab = "transitivity", xlim=c(0.25, 0.5)); abline(v = transitivity(g), col = "red")
hist(unlist(null.distr[,2]), xlab = "in-degree"); abline(v = sd(degree(g, mode = "in")), col = "red")
hist(unlist(null.distr[,3]), xlab = "out-degree", xlim=c(0, 22)); abline(v = sd(degree(g, mode = "out")), col = "red")
hist(unlist(null.distr[,4]), xlab = "density"); abline(v = edge_density(g), col = "red")

mean(transitivity(g) > null.distr[, 1])
mean(sd(degree(g, mode = "in")) > null.distr[, 2])
mean(sd(degree(g, mode = "out")) > null.distr[, 3])
mean(edge_density(g) > null.distr[, 4])

# simuliamo anche mod6_gwdsp per capire quali sono le differenze tra mod6_gwdsp e mod6_gwdsp_2
sim = simulate(mod6_gwdsp, nsim = 100, verbose = TRUE, seed = 1)

null.distr = matrix(,100,4)
for(b in 1:100){
  null.distr[b,]  = fnc(sim[[b]])
}
dev.new()
par(mfrow = c(4,1))
hist(unlist(null.distr[,1]), xlab = "transitivity"); abline(v = transitivity(g), col = "red")
hist(unlist(null.distr[,2]), xlab = "in-degree", ); abline(v = sd(degree(g, mode = "in")), col = "red")
hist(unlist(null.distr[,3]), xlab = "out-degree"); abline(v = sd(degree(g, mode = "out")), col = "red")
hist(unlist(null.distr[,4]), xlab = "density"); abline(v = edge_density(g), col = "red")

mean(transitivity(g) > null.distr[, 1])
mean(sd(degree(g, mode = "in")) > null.distr[, 2])
mean(sd(degree(g, mode = "out")) > null.distr[, 3])
mean(edge_density(g) > null.distr[, 4])

# simuliamo anche mod0 pe divertimento
sim0 = simulate(mod0, nsim = 100, verbose = TRUE, seed = 1)

null.distr = matrix(,100,4)
for(b in 1:100){
  null.distr[b,]  = fnc(sim0[[b]])
}
dev.new()
par(mfrow = c(4,1))
hist(unlist(null.distr[,1]), xlab = "transitivity"); abline(v = transitivity(g), col = "red")
hist(unlist(null.distr[,2]), xlab = "in-degree", ); abline(v = sd(degree(g, mode = "in")), col = "red")
hist(unlist(null.distr[,3]), xlab = "out-degree"); abline(v = sd(degree(g, mode = "out")), col = "red")
hist(unlist(null.distr[,4]), xlab = "density"); abline(v = edge_density(g), col = "red")

mean(transitivity(g) > null.distr[, 1])
mean(sd(degree(g, mode = "in")) > null.distr[, 2])
mean(sd(degree(g, mode = "out")) > null.distr[, 3])
mean(edge_density(g) > null.distr[, 4])

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

prMaxes = c(80)
for (i in 1:80) {
  prMaxes[i] <- max(sbm1$probMemberships[i,])
}
order(prMaxes, decreasing = FALSE)
prMaxes[order(prMaxes, decreasing = FALSE)]

# Prendiamo solo quelli sotto il 0.95
ordSBM1 = order(prMaxes, decreasing = FALSE)[1:7]
prMaxes[ordSBM1]
m = matrix(, nrow = 7, ncol = 6, dimnames = list(V(g)$name[ordSBM1], c("1", "2", "3", "4", "5", "6")))
for (i in 1:7) {
  m[i,1] = sbm1$probMemberships[ordSBM1[i], 1]
  m[i,2] = sbm1$probMemberships[ordSBM1[i], 2]
  m[i,3] = sbm1$probMemberships[ordSBM1[i], 3]
  m[i,4] = sbm1$probMemberships[ordSBM1[i], 4]
  m[i,5] = sbm1$probMemberships[ordSBM1[i], 5]
  m[i,6] = sbm1$probMemberships[ordSBM1[i], 6]
}
m

membership_colors = palette.colors(9)[-1]
plot(g, layout = l, rescale = FALSE, vertex.size=10, vertex.label.font = 2, vertex.color = membership_colors[sbm1$memberships])
legend("topright", legend = c("1", "2", "3", "4", "5", "6"), fill = membership_colors)
