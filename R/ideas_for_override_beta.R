beta <- matrix(0, nrow = 2, ncol = 4)
rownames(beta) <- c("L0", "A0")
colnames(beta) <- paste0("N", 0:3)

override_beta <- list("N2>=1" = c("N1" = 1.2, "N3" = 0.8),"N2>=2" = c("N1" = 0.6, "N3" = -0.2),
                      "A0" = c("N1" = 1, "N0" = 3, "N3" = -1))

for (bb in 1:length(override_beta)) {
  if (names(override_beta)[bb] %in% rownames(beta)) {
    beta[names(override_beta)[bb], names(override_beta[[bb]])] <- override_beta[[bb]]
  } else {
    beta <- rbind(beta, matrix(0, nrow = 1, ncol = ncol(beta)))
    beta[nrow(beta), names(override_beta[[bb]])] <- override_beta[[bb]]
    rownames(beta)[nrow(beta)] <- names(override_beta)[bb]
  }
}

L0 <- 0.6
A0 <- 0.5
N0 <- 0
N1 <- 0
N2 <- 2
N3 <- 0

mat <- do.call("rbind", lapply(1:nrow(beta), function(bb) {
  eval(parse(text = rownames(beta)[bb])) %*% beta[bb,]
}))

colSums(mat)
eval(parse(text = rownames(beta)[2]))


phi <- function(i) {
  do.call("rbind", lapply(1:nrow(beta), function(bb) {
    rowname_expr <- rownames(beta)[bb]

    # Evaluate the expression using the i-th element
    value <- eval(parse(text = rowname_expr))[i]

    # Perform matrix multiplication
    value %*% beta[bb,]
  }))
}
