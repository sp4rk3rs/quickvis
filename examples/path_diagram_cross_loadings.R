#### Fit reproducible population models illustrating cross-loadings and a bifactor
library(quickvis)
library(lavaan)

#### x3 is measured by both F1 and F2
lambda <- cbind(F1 = c(.8, .7, .6, 0, 0, 0), F2 = c(0, 0, .3, .8, .7, .6))
rownames(lambda) <- paste0("x", 1:6)
phi <- matrix(c(1, .3, .3, 1), 2)
sigma <- lambda %*% phi %*% t(lambda) + diag(rep(.5, 6))
cross_fit <- cfa('F1 =~ x1 + x2 + x3
                 F2 =~ x4 + x5 + x6 + x3',
                sample.cov = sigma, sample.nobs = 1000,
                sample.cov.rescale = FALSE, std.lv = TRUE)
path_diagram(cross_fit)
path_diagram(cross_fit, variance = "x3", tracing = TRUE)

#### General factor g plus three specific factors, each with three indicators
lambda <- cbind(g = rep(.7, 9), S1 = c(.5, .4, .6, rep(0, 6)),
                S2 = c(rep(0, 3), .5, .4, .6, rep(0, 3)),
                S3 = c(rep(0, 6), .5, .4, .6))
rownames(lambda) <- paste0("x", 1:9)
sigma <- tcrossprod(lambda) + diag(rep(.5, 9))
bifactor_fit <- cfa('g =~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9
                    S1 =~ x1 + x2 + x3
                    S2 =~ x4 + x5 + x6
                    S3 =~ x7 + x8 + x9',
                   sample.cov = sigma, sample.nobs = 1000,
                   sample.cov.rescale = FALSE, std.lv = TRUE, orthogonal = TRUE)
path_diagram(bifactor_fit) # Automatically places g below the indicators
path_diagram(bifactor_fit, layout = "bifactor", general = "g")
path_diagram(bifactor_fit, covariance = c("x1", "x2"), tracing = TRUE)
path_diagram(bifactor_fit, layout = "tree") # All factors above the indicators
