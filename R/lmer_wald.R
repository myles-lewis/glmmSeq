
lmer_wald <- function(fixef, hyp.matrix.1, hyp.matrix.2, vcov.) {
  
  hyp.list <- lapply(seq_along(hyp.matrix.1), function(i) {
    # source car:::Anova.II.mer
    hyp.matrix.term <- if (nrow(hyp.matrix.1[[i]]) == 0) 
    {hyp.matrix.2[[i]]
    } else t(car_ConjComp(t(hyp.matrix.1[[i]]), t(hyp.matrix.2[[i]]), vcov.))
    hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 
                                              1, function(x) all(x == 0)), , drop = FALSE]
    hyp.matrix.term
  })
  # source car:::linearHypothesis.mer()
  b <- fixef
  V <- vcov.
  # rhs = 0
  chi_val <- lapply(hyp.list, function(L) {
    as.vector(t(L %*% b) %*% solve(L %*% V %*% 
                                     t(L)) %*% (L %*% b))
  })
  df <- unlist(lapply(hyp.list, NROW))
  list(chisq = setNames(unlist(chi_val), names(hyp.matrix.1)),
       df = df)
}



## Source car:::ConjComp
car_ConjComp <- function (X, Z = diag(nrow(X)), ip = diag(nrow(X))) 
{
  xq <- qr(t(Z) %*% ip %*% X)
  if (xq$rank == 0) 
    return(Z)
  Z %*% qr.Q(xq, complete = TRUE)[, -(1:xq$rank)]
}


## Source car:::relatives
car_relatives <- function (term, names, factors) 
{
  is.relative <- function(term1, term2) {
    all(!(factors[, term1] & (!factors[, term2])))
  }
  if (length(names) == 1) 
    return(NULL)
  which.term <- which(term == names)
  (1:length(names))[-which.term][sapply(names[-which.term], 
                                        function(term2) is.relative(term, term2))]
}
