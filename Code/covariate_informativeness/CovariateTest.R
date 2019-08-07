require(splines)
require(matrixStats)
require(DescTools)

CovTest.n <- function (pvalue, covariate, cutoffs = quantile(pvalue, c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2)),
		grps = c(2, 4, 8, 16, 32), perm.no = 999, n.max = 100000, silence = TRUE) {
	
	save.seed <- get(".Random.seed", .GlobalEnv)
	
	n <- length(pvalue)
	if (n > n.max) {
		ind <- sample(1:n, n.max)
		pvalue <- pvalue[ind]
		covariate <- covariate[ind]
	}
	
	
	stat.p.a <- stat.p.a1 <- stat.p.a2 <- array(NA, c(length(cutoffs), length(grps), perm.no))
	stat.o.a <- stat.o.a1 <- stat.o.a2 <-array(NA, c(length(cutoffs), length(grps)))
	
	for (i in 1:length(cutoffs)) {
		cutoff <- cutoffs[i]
		x <- as.numeric(pvalue <= cutoff) 
		for (j in 1:length(grps)) {
			if (!silence) cat('.')
			grp <- grps[j]
			y <- cut(covariate, c(min(covariate) - 0.01, quantile(covariate, 1 / grp * (1:(grp-1))), max(covariate) + 0.01))
			
			if (!silence) cat('.')
			mat <- table(x, y)
			test.obj1 <- chisq.test(mat)
			test.obj2 <- CochranArmitageTest(mat)
			stat.o.a1[i, j] <- -pchisq(test.obj1$statistic, df = test.obj1$parameter, lower.tail = FALSE, log.p = TRUE)
			stat.o.a2[i, j] <- -log(2) - pnorm(abs(test.obj2$statistic), lower = FALSE, log.p = TRUE)
			stat.o.a[i, j] <- max(c(stat.o.a1[i, j], stat.o.a2[i, j]))
			
			assign(".Random.seed", save.seed, .GlobalEnv)
			
			lpvs <- sapply(1:perm.no, function (k) {
						xp <- sample(x)		
						mat <- table(xp, y)
						test.obj1 <- chisq.test(mat)
						test.obj2 <- CochranArmitageTest(mat)
						c(-pchisq(test.obj1$statistic, df = test.obj1$parameter, lower.tail = FALSE, log.p = TRUE),
								-log(2) - pnorm(abs(test.obj2$statistic), lower = FALSE, log.p = TRUE))
					})
			
			stat.p.a1[i, j, ]  <- lpvs[1, ]
			stat.p.a2[i, j, ]  <- lpvs[2, ]
			stat.p.a[i, j, ] <- pmax(lpvs[1, ], lpvs[2, ])
		}
	}
	
	stat.o1 <- colMaxs(stat.o.a)
	stat.o2 <- rowMaxs(stat.o.a)
	stat.o <- max(stat.o1)
	x.cut.optim <- grps[which.max(stat.o1)]
	p.cut.optim <- cutoffs[which.max(stat.o2)]
	stat.p <- apply(stat.p.a, 3, max)
	p.value <- mean(c(stat.p, stat.o) >= stat.o)
	
	stat.o1 <- max(stat.o.a1)
	stat.p1 <- apply(stat.p.a1, 3, max)
	p.value1 <- mean(c(stat.p1, stat.o1) >= stat.o1)
	
	stat.o2 <- max(stat.o.a2)
	stat.p2 <- apply(stat.p.a2, 3, max)
	p.value2 <- mean(c(stat.p2, stat.o2) >= stat.o2)
	
	return(list(stat.o = stat.o,  stat.p = stat.p, p.value = p.value,
					stat.o1 = stat.o1,  stat.p1 = stat.p1, p.value1 = p.value1, 
					stat.o2 = stat.o2,  stat.p2 = stat.p2, p.value2 = p.value2,
					x.cut.optim = x.cut.optim,  p.cut.optim = p.cut.optim))
}



CovTest.c <- function (pvalue, covariate, cutoffs = quantile(pvalue, c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2)), 
		perm.no = 999, n.max = 100000, silence = TRUE) {
	
	save.seed <- get(".Random.seed", .GlobalEnv)
	
	n <- length(pvalue)
	if (n > n.max) {
		ind <- sample(1:n, n.max)
		pvalue <- pvalue[ind]
		covariate <- covariate[ind]
	}
	
	y <- factor(covariate)
	
	stat.p.a <- array(NA, c(length(cutoffs), perm.no))
	stat.o.a <- array(NA, c(length(cutoffs)))
	
	for (i in 1:length(cutoffs)) {
		
		assign(".Random.seed", save.seed, .GlobalEnv)
		
		cutoff <- cutoffs[i]
		x <- as.numeric(pvalue <= cutoff) 
		
		if (!silence) cat('.')
		test.obj <- chisq.test(table(x, y))
		stat.o.a[i] <- -pchisq(test.obj$statistic, df = test.obj$parameter, lower.tail = FALSE, log.p = TRUE)
		stat.p.a[i, ] <- sapply(1:perm.no, function (k) {
					xp <- sample(x)		
					test.obj.p <- chisq.test(table(xp, y))
					-pchisq(test.obj.p$statistic, df = test.obj.p$parameter, lower.tail = FALSE, log.p = TRUE)
				})
		
	}
	
	
	stat.o <- max(stat.o.a)
	p.cut.optim <- cutoffs[which.max(stat.o.a)]
	stat.p <- apply(stat.p.a, 2, max)
	p.value <- mean(c(stat.p, stat.o) >= stat.o)
	return(list(stat.o.a = stat.o.a, stat.o = stat.o, stat.p = stat.p, p.value = p.value, p.cut.optim = p.cut.optim))
}

## Example 1
#n <- 450000
#covariate <- rnorm(n)
#eta <- 3.5 + scale(ns(covariate, df = 4) %*% rnorm(4)) * 0.25
#pi <- exp(eta) / (1 + exp(eta))
#mean(pi)
#
#truth <- rbinom(n, 1 - pi, size = 1)
#
#pvalue <- pnorm(ifelse(truth == 1, rnorm(n, mean = 2), rnorm(n)), lower = FALSE)
#boxplot(-log(pvalue) ~ truth)
#hist(pvalue)
#
#obj <- CovTest.n(pvalue, covariate, perm.no = 999)
#obj$p.value
#obj$stat.o
#
#
## Example 2
#n <- 10000
#covariate <- gl(5, n / 5)
#eta <- 3.5 + scale(model.matrix(~ covariate)[, -1, drop = FALSE] %*% rnorm(4)) * 0.25
#pi <- exp(eta) / (1 + exp(eta))
#mean(pi)
#
#truth <- rbinom(n, 1 - pi, size = 1)
#
#pvalue <- pnorm(ifelse(truth == 1, rnorm(n, mean = 2), rnorm(n)), lower = FALSE)
#boxplot(-log(pvalue) ~ truth)
#hist(pvalue)
#
#obj <- CovTest.c(pvalue, covariate, perm.no = 999)
#obj$stat.o
#
