###############################################
###############################################
## script from Derkarabetian S., Castillo S., Peter K.K., Ovchinnikov S., Hedin M. "An Empirical Demonstration of Unsupervised Machine Learning in Species Delimitation"
###############################################
###############################################


#required packages
library("adegenet")
library("randomForest")
library("PCDimension")
library("mclust")
library("cluster")
library("MASS")
library("factoextra")
library("tsne")

## This assumes your input file is structure/adegenet format (.str). You can also import data in .csv format, although
## depending on data type, you will have to convert:
## If importing .csv with raw nucleotides/haplotypes, convert to factor, go to random forest step.
## If importing .csv with nucleotide data in one-hot format, covert to numeric, go to PCA step.

###############################################
###############################################
# PCA and DAPC
###############################################
###############################################

# import str file. Adjust input file name, n.ind, and n.loc for specific file/dataset.
# example dataset used in this study
# data <- read.structure("Metano_UCE_SNPs_70percent_random_structure-adegenet.str", n.ind=30, n.loc=316, onerowperind=FALSE, col.lab=1, col.pop=3, col.others=NULL, row.marknames=NULL, NA.char="-9", pop=NULL, ask=FALSE, quiet=FALSE)

# data <- read.structure("input.str", n.ind=XX, n.loc=XX, onerowperind=FALSE, col.lab=1, col.pop=0, col.others=NULL, row.marknames=NULL, NA.char="-9", pop=NULL, ask=FALSE, quiet=FALSE)
data_scaled <- scaleGen(data, center=FALSE, scale=FALSE, NA.method=c("zero"), nf)

# PCA, can adjust nf to include more components
pca1 <- dudi.pca(data_scaled, center=TRUE, scale=TRUE, scannf=FALSE, nf=2)
plot(pca1$li, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="PCA", pch=16)
s.label(pca1$li, clabel=0.5, grid=0)

# DAPC (interactive, requires input)
# max.n.clust equal to number of pops
clusters <- find.clusters(data, max.n.clust=10, n.iter=1e6, n.start=10)
results <- dapc(data, grp=clusters$grp, perc.pca=NULL)
assignplot(results)
compoplot(results)
grp_k <- nlevels(clusters$grp)

# PCA with DAPC groups
plot(pca1$li, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="PCA with DAPC clusters", col=results$grp, pch=16)



###############################################
###############################################
# into the Random Forest, unsupervised
###############################################
###############################################

# convert genind scaled data to factors for randomForest
data_conv <- as.data.frame(data_scaled)
data_conv[is.na(data_conv)] <- ""
data_conv[sapply(data_conv, is.integer)] <- lapply(data_conv[sapply(data_conv, is.integer)], as.factor)
data_conv[sapply(data_conv, is.character)] <- lapply(data_conv[sapply(data_conv, is.character)], as.factor)
nsamp <- nrow(data_conv)

# unsupervised random forest
rftest <- randomForest(data_conv, ntree=5000)

###############
# classic MDS
###############

# cMDS with optimal number of components to retain using broken-stick
# may need to adjust number of dimensions if given error
cmdsplot1 <- MDSplot(rftest, results$grp, nsamp-1)
cmdsplot_bstick <- bsDimension(cmdsplot1$eig)
cmdsplot2 <- MDSplot(rftest, results$grp, cmdsplot_bstick)

# cMDS with optimal DAPC k and clusters
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="cMDS DAPC optimal K and clusters", col=results$grp, pch=16)

# pam clustering on proximity scores with optimal k from DAPC
DAPC_pam_clust_prox <- pam(rftest$proximity, grp_k)
# cMDS with optimal k of DAPC and clusters via PAM
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="cMDS DAPC optimal K and clusters (PAM clustering)", col=DAPC_pam_clust_prox$clustering, pch=16)
# mean silhouette width of k
# can test multiple k values and select optimal k with highest value
mean(silhouette(DAPC_pam_clust_prox)[, "sil_width"])

# pam clustering on cMDS output with optimal k from DAPC
DAPC_pam_clust_cMDS <- pam(cmdsplot1$points, grp_k)
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="cMDS DAPC optimal K and clusters (PAM clustering)", col=DAPC_pam_clust_cMDS$clustering, pch=16)
# mean silhouette width of k
# can test multiple k values and select optimal k with highest value
mean(silhouette(DAPC_pam_clust_cMDS)[, "sil_width"])

# determine optimal k from cMDS using gap statistic with PAM clusters from proximity scores
# can adjust k.max
cmds_nbclust <- fviz_nbclust(cmdsplot1$points, kmeans, nstart = 25,  method = "gap_stat", nboot = 500) + labs(subtitle = "Gap statistic method")
cmds_nbclust
cmds_nbclust_k <- cmds_nbclust[["layers"]][[4]][["data"]][["xintercept"]]
# pam clustering with optimal k from gap statistic
cmds_nbclust_clust <- pam(cmdsplot1$points, cmds_nbclust_k)
# cMDS with optimal k of RF via gap statistic and clusters via PAM (euc)
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="proximity RF gap statistic optimal K and clusters (PAM clustering)", col=cmds_nbclust_clust$clustering, pch=16)

# determine optimal k from cMDS via hierarchical clustering with BIC
# adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
cmdsplot_clust <- Mclust(cmdsplot2$points)
mclust_grps_cmdsplot <- as.numeric(cmdsplot_clust$classification)
max(mclust_grps_cmdsplot)
# cMDS with optimal k and clusters of RF via hierarchical clustering
plot(cmdsplot2$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="cMDS RF optimal K and clusters (hierarchical clustering)", col=mclust_grps_cmdsplot, pch=16)
mclust_grps_cmdsplot

s.label(cmdsplot2$points, clabel=0.5, grid=0)

###############
# isotonic MDS
###############

# isoMDS
isomdsplot <- isoMDS(1-rftest$proximity)
# "The output of cmdscale on 1 - rf$proximity is returned invisibly" (MDSplot documentation)
plot(isomdsplot$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="isoMDS DAPC optimal K and clusters", col=results$grp, pch=16)

# pam clustering with optimal k from DAPC
DAPC_pam_clust_iso <- pam(isomdsplot$points, grp_k)
# cMDS with optimal k of DAPC and clusters via PAM
plot(isomdsplot$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="cMDS DAPC optimal K and clusters (PAM clustering)", col=DAPC_pam_clust_iso$clustering, pch=16)
# mean silhouette width of k
# can test multiple k values and select k with hiest value
mean(silhouette(DAPC_pam_clust_iso)[, "sil_width"])

# determine optimal k using gap statistic
# can adjust k.max
isomds_nbclust <- fviz_nbclust(isomdsplot$points, kmeans, nstart = 25,  method = "gap_stat", nboot = 500) + labs(subtitle = "Gap statistic method")
isomds_nbclust
isomds_nbclust_k <- isomds_nbclust[["layers"]][[4]][["data"]][["xintercept"]]
# pam clustering with optimal k from gap statistic
isomds_nbclust_clust2 <- pam(rftest$proximity, isomds_nbclust_k)
# isoMDS with optimal k of RF via gap statistic and clusters via PAM
plot(isomdsplot$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="isoMDS RF gap statistic optimal K and clusters (PAM clustering)", col=isomds_nbclust_clust2$clustering, pch=16)

# determine optimal k of RF via hierarchical clustering with BIC
# adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
isomdsplot_clust <- Mclust(isomdsplot$points)
mclust_grps_isomdsplot2 <- as.numeric(isomdsplot_clust$classification)
max(mclust_grps_isomdsplot2)
# isoMDS with optimal k and clusters of RF via hierarchical clustering
plot(isomdsplot$points, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="isoMDS RF optimal K and clusters (hierarchical clustering)", col=mclust_grps_isomdsplot2, pch=16)
mclust_grps_isomdsplot2

s.label(isomdsplot$points, clabel=0.5, grid=0)



###############################################
###############################################
# t-SNE
###############################################
###############################################

# prepare plot labels and such
# this assumes you have a population assignment (a priori species) column in the .str file.
colors = rainbow(length(unique(data$pop)))
names(colors) = unique(data$pop)
ecb = function(x,y){plot(x,t='n'); text(x, labels=data$pop, col=colors[data$pop])}
# OR
# this makes it so it is grouped by DAPC clusters
colors = rainbow(length(unique(results$grp)))
names(colors) = unique(results$grp)
ecb = function(x,y){plot(x,t='n'); text(x, labels=results$grp, col=colors[results$grp])}

# t-SNE on principal components of scaled data
# adjust perplexity, initial_dims
# can do k=3 for 3D plot
# should do only <50 variables
# can do it on pca$li (if you reduce the number of components), or on cmdsplot2$points
tsne_p5 = tsne(pca1$tab, epoch_callback=ecb, max_iter=5000, perplexity=5, initial_dims=5)

# tSNE plot with DAPC groups
plot(tsne_p5, main="t-SNE perplexity=5 with DAPC optimal k and clusters", col=results$grp, pch=16)
s.label(tsne_p5, clabel=0.5, grid=0)

# pam clustering with optimal k from DAPC
DAPC_pam_clust2_tsne <- pam(tsne_p5, grp_k)
# mean silhouette width of k
# can test multiple k values and select k with highest value
mean(silhouette(DAPC_pam_clust2_tsne)[, "sil_width"])

# clustering for perplexity=5
# determine optimal k using gap statistic
# can adjust k.max
tsne_p5_nbclust <- fviz_nbclust(tsne_p5, kmeans, nstart = 25,  method = "gap_stat", nboot = 500) + labs(subtitle = "Gap statistic method")
tsne_p5_nbclust
tsne_p5_nbclust_k <- tsne_p5_nbclust[["layers"]][[4]][["data"]][["xintercept"]]
# pam clustering with optimal k from gap statistic
tsne_p5_nbclust_clust <- pam(tsne_p5, tsne_p5_nbclust_k)
# t-SNE with optimal k of RF via gap statistic and clusters via PAM
plot(tsne_p5, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="t-SNE p5 RF gap statistic optimal K and clusters (PAM clustering)", col=tsne_p5_nbclust_clust$clustering, pch=16)

# determine optimal k of RF via hierarchical clustering with BIC
# adjust G option to reasonable potential cluster values, e.g. for up to 12 clusters, G=1:12
tsne_p5_clust <- Mclust(tsne_p5)
mclust_grps_tsne_p5 <- as.numeric(tsne_p5_clust$classification)
max(mclust_grps_tsne_p5)
# t-SNE p5 with optimal k and clusters of RF via hierarchical clustering
plot(tsne_p5, xlab="Scaling Dimension 1", ylab="Scaling Dimension 2", main="t-SNE p5 RF optimal K and clusters (hierarchical clustering)", col=mclust_grps_tsne_p5, pch=16)
mclust_grps_tsne_p5


