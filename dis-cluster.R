##################################################
# Discovery cohort                               #
# Sparse multiple canonical correlation analysis #
# + K-means + NMF meta-clustering                #
##################################################
library(PMA)
library(NMF)
set.seed(1234)

##################################################
# Toy data: three brain regions                  #
##################################################
n_samples     <- 459
n_genes_DLPFC <- 18629
n_genes_AC    <- 19147
n_genes_PCC   <- 19017

# Latent subject factor
u <- matrix(rnorm(n_samples), ncol = 1)

# Region-specific loading patterns
v_DLPFC <- matrix(c(rep(1,  25), rep(0, n_genes_DLPFC - 25)), ncol = 1)
v_AC    <- matrix(c(rep(0.5,25), rep(0, n_genes_AC    - 25)), ncol = 1)
v_PCC   <- matrix(c(rep(0.5,25), rep(0, n_genes_PCC   - 25)), ncol = 1)

# Simulated expression matrices (samples x genes)
gx_DLPFC <- u %*% t(v_DLPFC) + matrix(rnorm(n_samples * n_genes_DLPFC), nrow = n_samples)
gx_AC    <- u %*% t(v_AC)    + matrix(rnorm(n_samples * n_genes_AC),    nrow = n_samples)
gx_PCC   <- u %*% t(v_PCC)   + matrix(rnorm(n_samples * n_genes_PCC),   nrow = n_samples)

##################################################
# Sparse multiple canonical correlation analysis #
##################################################
xlist <- list(gx_DLPFC, gx_AC, gx_PCC)

# Choose penalties via permutation (low nperms here for illustration)
perm.out <- MultiCCA.permute(xlist, nperms = 10, type = "standard")

cca_fit <- MultiCCA(xlist, type="standard", penalty=perm.out$bestpenalties, ws=perm.out$ws.init, ncomponents=10, standardize = TRUE)

# Canonical weight vectors
w_DLPFC <- cca_fit$ws[[1]]
w_AC    <- cca_fit$ws[[2]]
w_PCC   <- cca_fit$ws[[3]]

# Canonical component scores (samples x components)
cv_DLPFC <- gx_DLPFC %*% w_DLPFC
cv_AC    <- gx_AC    %*% w_AC
cv_PCC   <- gx_PCC   %*% w_PCC

#########################################
# K-means clustering within each region #
#########################################
region_data <- list(
  DLPFC = cv_DLPFC,
  AC    = cv_AC,
  PCC   = cv_PCC
)

clusters  <- list()
centroids <- list()
for (region in names(region_data)) {
  x    <- region_data[[region]]
  kfit <- kmeans(x, centers = 2)
  clusters[[region]]  <- data.frame(cluster = kfit$cluster)
  centroids[[region]] <- kfit$centers
}

##################################################
# Meta-clustering across regions using NMF       #
##################################################
# Binary cluster membership matrix:
#   rows = region-specific clusters
#   cols = subjects (aligned by rownames)
X <- rbind(
  DLPFC1 = as.numeric(clusters$DLPFC$cluster == 1),
  DLPFC2 = as.numeric(clusters$DLPFC$cluster == 2),
  AC1    = as.numeric(clusters$AC$cluster    == 1),
  AC2    = as.numeric(clusters$AC$cluster    == 2),
  PCC1   = as.numeric(clusters$PCC$cluster   == 1),
  PCC2   = as.numeric(clusters$PCC$cluster   == 2)
)
colnames(X) <- rownames(clusters$DLPFC)

# NMF to identify shared “meta-clusters” across regions
k_meta  <- 2
nmf_fit <- nmf(X, rank = k_meta, method = "lee", seed = 123456, nrun = 10)

# Meta-cluster assignment per subject
H <- coef(nmf_fit)                      # meta-clusters x subjects
meta_cluster <- apply(H, 2, which.max)  # length = n_samples
names(meta_cluster) <- colnames(H)
