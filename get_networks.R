# Analysis pipeline extracting sparse, robust brain subnetworks predictive of traits.
# Steven Winter, Sept. 2023
#
# 1. Select optimal sparsity parameters with cross validation.
# 2. Fit optimal model on replicates of brain data.
# 3. Compute theoretical false discovery bounds.
# Note: steps 1-2 are computationally intensive and may require cluster computing. 
# This code is modular to allow easy parallelization.
#
# References:
# Variable selection with error control: Another look at Stability Selection, Rajen D. Shah, Richard J. Samworth
# Network classification with applications to brain connectomics, Arroyo Relión et al 
set.seed(1234)
library(dplyr)
library(tidyverse)

library(parallel)
library(devtools)
#install_github("jesusdaniel/graphclass")
library(graphclass)

# Connectomes should be a real valued array with dimensions (V, V, n)
# response should be a binary vector of length n

setwd('/Users/alex/AlexBadea_MyCodes/CVN_GraphClass/');
input_path <- paste0(getwd(), '/'); 


cv_path <- paste0(input_path, 'cv/')
stability_path<- paste0(input_path , "/stability/")

if (!dir.exists(cv_path)){
  dir.create(cv_path)
} else {
  print("Dir already exists!")
}

if (!dir.exists(stability_path)){
  dir.create(stability_path)
} else {
  print("Dir already exists!")
}



#connectomes <- readRDS(paste0(input_path, "connectomes.rds"))

#response
load(paste0(input_path, "/response.rda"))
load(paste0(input_path, "/connectivity_plain.rda"))
dim(connectivity)
dim(response)

index = which(response$treatment %in% c("treadmill","wheel_only"))

#ind = ! response$treatment=="sedentary"

Y= response[index,]$treatment
Y[Y=="treadmill"] =1 ; Y[Y=="wheel_only"] =-1
Y = as.numeric(Y)  
connectomes =  connectivity[,,index]
dim(connectomes)  

                 

V <- dim(connectomes)[1]
n <- dim(connectomes)[3]
p <- V*(V-1)/2
X <- array(dim=c(n,p))
for(i in 1:n){
  X[i,] <- matrix_to_vec(connectomes[,,i])
}

#################################################################################################
####################################### CROSS VALIDATION ########################################
#################################################################################################

fold_index <- (1:length(Y) %% 10) + 1
fold_index <- (1:length(Y) %% 4) + 1

fold_index <- sample(fold_index, replace = FALSE)

# library(caret)
# trainIndex <- createDataPartition(Y, p = .8)


DV <- construct_D(nodes = V) 
gclist <- list()
lams <- 10^seq(-4,2,by=0.5)
rhos <- 10^seq(-2,2,by=0.5)

#lams <- 10^seq(-4,-3.5,by=0.5)
#rhos <- 10^seq(-2,-1.5,by=0.5)


# lams <- 10^seq(-4,-2,by=1)
# rhos <- 10^seq(-2,-2,by=1)
# #

# #start unparalelled cross validation
# for(lambda in lams){
#   for(rho in rhos){
#     #for(fold in 1:10){
#     for(fold in 1:2){
# 
#       print(fold)
#       foldout <- which(fold_index == fold)
#       gclist[[fold]] <- graphclass(X = X[-foldout,], Y = factor(Y[-foldout]),
#                                   Xtest = X[foldout,], Ytest = factor(Y[foldout]),
#                                   type = "intersection",
#                                   lambda = lambda, rho = rho, gamma = 1e-5,
#                                   D = DV)
#       }
#     errors <- unlist(lapply(gclist, function(gc) gc$test_error))
#     gc_full <- graphclass(X = X, Y = factor(Y),
#                          lambda = lambda, rho = rho, gamma = 1e-5,
#                          D = DV)
#     full_beta <- gc_full$beta
#     num_edges <- sum(gc_full$beta != 0)
#     # Save the hyperparameters, the average error, standard deviation,
#     # number of selected edges, and full coefficients
#     results <- list(lambda=lambda,
#                    rho=rho,
#                    mean_error = mean(errors),
#                    error_sd = sd(errors),
#                    num_edges = num_edges,
#                    full_beta = full_beta)
#     fname <- paste0("cv_",lambda,"_",rho, ".rds")
#     saveRDS(results, file=paste0(cv_path, fname))
#   }
# }
# 
# #end unparalelled cross validation




#start paralell cross validation
print('start paralell')
mclapply(1:length(lams), function(i) {
#mclapply(1:1, function(i) {
    
  lambda=lams[i]

#for(lambda in lams){
  for(rho in rhos){
    #for(fold in 1:10){
    for(fold in 1:4){
      
      print(fold)
      foldout <- which(fold_index == fold)
      gclist[[fold]] <- graphclass(X = X[-foldout,], Y = factor(Y[-foldout]),
                                   Xtest = X[foldout,], Ytest = factor(Y[foldout]),
                                   type = "intersection",
                                   lambda = lambda, rho = rho, gamma = 1e-5,
                                   D = DV)
    }
    errors <- unlist(lapply(gclist, function(gc) gc$test_error))
    gc_full <- graphclass(X = X, Y = factor(Y),
                          lambda = lambda, rho = rho, gamma = 1e-5,
                          D = DV)
    full_beta <- gc_full$beta
    num_edges <- sum(gc_full$beta != 0)
    # Save the hyperparameters, the average error, standard deviation,
    # number of selected edges, and full coefficients
    results <- list(lambda=lambda,
                    rho=rho,
                    mean_error = mean(errors),
                    error_sd = sd(errors),
                    num_edges = num_edges,
                    full_beta = full_beta)
    fname <- paste0("cv_",lambda,"_",rho, ".rds")
    saveRDS(results, file=paste0(cv_path, fname))
  }

}  , mc.cores = 6)



#end paralell cross validation




# Load cross-validation results
# Sort by accuracy, break ties by choosing the sparsest model
fnames <- paste0(cv_path, list.files(path = cv_path))
raw_results <- lapply(fnames, readRDS)
results <- data.frame(lambdas = unlist(lapply(raw_results, function(x) x$lambda)),
                     rhos = unlist(lapply(raw_results, function(x) x$rho)),
                     means = unlist(lapply(raw_results, function(x) x$mean)),
                     sds = unlist(lapply(raw_results, function(x) x$error_sd)),
                     num_edges = unlist(lapply(raw_results, function(x) x$num_edges)),
                     fname = fnames)
results_sorted <- results %>% 
  arrange(means, num_edges)


# Save the best model
print(head(results_sorted))
best_file <- readRDS(results_sorted$fname[1])
saveRDS(best_file, file=paste0(input_path, "best.rds"))


# Save the second best model
#print(head(results_sorted))
#best_file <- readRDS(results_sorted$fname[8])
#saveRDS(best_file, file=paste0(input_path, "best.rds"))



#################################################################################################
###################################### STABILITY SELECTION ######################################
#################################################################################################
#stability_path <- "~/stability/"
#stability_path<-paste0(input_path, "stability/")
# Get optimal CV parameters for stability selection
best <- readRDS(paste0(input_path, "best.rds"))
lambda <- best$lambda
rho <- best$rho



B <- 50
#for(split in 1:B){
  
  mclapply(1:B, function(i) {
  split=i
    
  fold1 <- sort(sample(1:n, floor(n/2)))
  fold2 <- (1:n)[!(1:n %in% fold1)]
  gc_f1 <- graphclass(X = X[fold1,], Y = factor(Y[fold1]), type = "intersection",
                     lambda = lambda, rho = rho, gamma = 1e-5)
  gc_f2 <- graphclass(X = X[fold2,], Y = factor(Y[fold2]), type = "intersection",
                     lambda = lambda, rho = rho, gamma = 1e-5)
  betas <- cbind(gc_f1$beta, gc_f2$beta)
  fname <- paste0("stability_",split,".rds")
  saveRDS(betas, file=paste0(stability_path, fname))
#}

  }  , mc.cores = 12)

####
  
  
  



# Concatenate coefficients from B runs into a [2*B]x[V*(V-1)/2] matrix
fnames <- paste0(stability_path, list.files(path = stability_path))
raw_results <- lapply(fnames, readRDS)
combined <- t(do.call(cbind, raw_results))

# B is the number of pairs that were run (e.g., B=50 means 100 total models fit)
# p is the number of predictors
# q_hat is the average number of nonzero coefficients
# probs is the probability of each coefficient being selected
B <- length(fnames)
p <- dim(combined)[2]
q_hat <- mean(rowSums(combined != 0))
probs <- colMeans(combined != 0)

# theta defines low selection probability, taken to be the average by default
# This needs to be less than 1/sqrt(3) for unimodal bound to be valid
theta <- q_hat/p
print("Valid theta")
print(theta < 1/sqrt(3))

# tau is the cutoff for a high selection probability edge
# This needs to be above 3/4 for the unimodal bound to be valid.
tau <- 0.75 #0.75
print("Valid tau")
print(tau > 3/4)

print("Number of high selection probability edges:")
print(sum(probs > tau))
print("Upper bound on number of low selection probability edges:")
C <- 4*(1-tau+1/(2*B))/(1+1/B)
print(ceiling(C*theta*q_hat))

# Save HSP adjacency matrix, as well as the top 50 edges
best <- readRDS(paste0(input_path, "best.rds"))
full_beta <- best$full_beta
full_beta[probs <= tau] <- 0
A_full <- get_matrix(full_beta)
saveRDS(A_full, file=paste0(input_path, "A_full.rds"))
A_50 <- A_full
sorted <- sort(abs(A_50), decreasing = T)
A_50[abs(A_50) <= sorted[2*51]] <- 0
saveRDS(A_50, file=paste0(input_path, "A_50.rds"))



########
saveRDS(results_sorted, file=paste0(input_path, "results_sorted.rds"))
HSP <- data.frame( q_hat = q_hat, theta = theta, tau = tau, myHSP=print(sum(probs > tau)), bound = ceiling(C*theta*q_hat))
saveRDS(HSP, file=paste0(input_path, "HSP.rds"))
print(HSP)
#######

