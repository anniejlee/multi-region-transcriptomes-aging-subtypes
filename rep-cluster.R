##################################################
# Replication cohort                             #
# Projection + K-means assignment + NMF          #
##################################################
library(PMA)
library(NMF)
set.seed(123)

##################################################
# Toy data: three brain regions                  #
##################################################
n_DLPFC_rep <- 633
n_AC_rep    <- 272
n_PCC_rep   <- 202

n_genes_DLPFC <- 18629
n_genes_AC    <- 19147
n_genes_PCC   <- 19017

# Latent subject factors for replication
u_DLPFC_rep <- matrix(rnorm(n_DLPFC_rep), ncol = 1)
u_AC_rep    <- matrix(rnorm(n_AC_rep),    ncol = 1)
u_PCC_rep   <- matrix(rnorm(n_PCC_rep),   ncol = 1)

# Simulated expression using discovery v-patterns
gx_DLPFC_rep <- u_DLPFC_rep %*% t(v_DLPFC) + matrix(rnorm(n_DLPFC_rep * n_genes_DLPFC), nrow = n_DLPFC_rep)
gx_AC_rep    <- u_AC_rep    %*% t(v_AC)    + matrix(rnorm(n_AC_rep * n_genes_AC),       nrow = n_AC_rep)
gx_PCC_rep   <- u_PCC_rep   %*% t(v_PCC)   +  matrix(rnorm(n_PCC_rep * n_genes_PCC),     nrow = n_PCC_rep)

##################################################
# Project replication samples into CCA space     #
##################################################
k <- ncol(w_DLPFC)  # number of canonical components 
cv_DLPFC_rep <- gx_DLPFC_rep %*% w_DLPFC
cv_AC_rep    <- gx_AC_rep    %*% w_AC
cv_PCC_rep   <- gx_PCC_rep   %*% w_PCC

##################################################
# Assign replication samples to discovery clusters
# using nearest centroids                        #
##################################################
assign_to_centroid <- function(x, centers) {
  dists <- sapply(seq_len(nrow(centers)), function(k) {
    rowSums((x - matrix(centers[k, ], nrow = nrow(x), ncol = ncol(x), byrow = TRUE))^2)
  })
  max.col(-dists)  # index of min distance
}

cluster_DLPFC_rep <- data.frame(cluster = assign_to_centroid(cv_DLPFC_rep, centroids$DLPFC), row.names = rownames(cv_DLPFC_rep))
cluster_AC_rep    <- data.frame(cluster = assign_to_centroid(cv_AC_rep,    centroids$AC), row.names = rownames(cv_AC_rep))
cluster_PCC_rep   <- data.frame(cluster = assign_to_centroid(cv_PCC_rep,   centroids$PCC), row.names = rownames(cv_PCC_rep))

##################################################
# Multiview NMF meta-clustering                  #
##################################################
# Collect all subject IDs across regions
ID_all    <- c(rownames(cluster_DLPFC_rep),
               rownames(cluster_AC_rep),
               rownames(cluster_PCC_rep))
ID_unique <- sort(unique(ID_all))
m         <- length(ID_unique)

# Binary cluster membership matrix:
#   rows = region-specific clusters
#   cols = subjects
X_rep <- matrix(0, nrow = 6, ncol = m)
rownames(X_rep) <- c("DLPFC1","DLPFC2","AC1","AC2","PCC1","PCC2")
colnames(X_rep) <- ID_unique

for (i in seq_len(m)) {
  id <- ID_unique[i]
  
  if (id %in% rownames(cluster_DLPFC_rep)) {
    if (cluster_DLPFC_rep[id, "cluster"] == 1) X_rep["DLPFC1", i] <- 1
    if (cluster_DLPFC_rep[id, "cluster"] == 2) X_rep["DLPFC2", i] <- 1
  }
  if (id %in% rownames(cluster_AC_rep)) {
    if (cluster_AC_rep[id, "cluster"] == 1) X_rep["AC1", i] <- 1
    if (cluster_AC_rep[id, "cluster"] == 2) X_rep["AC2", i] <- 1
  }
  if (id %in% rownames(cluster_PCC_rep)) {
    if (cluster_PCC_rep[id, "cluster"] == 1) X_rep["PCC1", i] <- 1
    if (cluster_PCC_rep[id, "cluster"] == 2) X_rep["PCC2", i] <- 1
  }
}

# NMF meta-clustering
k_meta <- 2
nmf_fit_rep <- nmf(X_rep, rank = k_meta, method = "lee", seed = 123456, nrun = 10)

# Meta-cluster assignment
H_rep <- coef(nmf_fit_rep)
meta_cluster_rep <- apply(H_rep, 2, which.max)
names(meta_cluster_rep) <- colnames(X_rep)
