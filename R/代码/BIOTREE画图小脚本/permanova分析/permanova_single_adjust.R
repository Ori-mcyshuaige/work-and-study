argv <- commandArgs(T)
if(length(argv) == 0){stop("Rscript permanova_single_adjust.R [input1] [input2] [prefix]
input1: Distance Matrix
input2: Phenotype
prefix: eg: Bray Euclid")
}


adonis1<-function (formula, data = NULL, permutations = 999, method = "bray", 
    strata = NULL, contr.unordered = "contr.sum", contr.ordered = "contr.poly") 
{

    TOL <- 1e-07
    Terms <- terms(formula, data = data)
    lhs <- formula[[2]]
    lhs <- eval(lhs, data, parent.frame())
    formula[[2]] <- NULL
    rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE)
    op.c <- options()$contrasts
    options(contrasts = c(contr.unordered, contr.ordered))
    rhs <- model.matrix(formula, rhs.frame)
    options(contrasts = op.c)
    grps <- attr(rhs, "assign")
    qrhs <- qr(rhs)
    rhs <- rhs[, qrhs$pivot, drop = FALSE]
    rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
    grps <- grps[qrhs$pivot][1:qrhs$rank]
    u.grps <- unique(grps)
    nterms <- length(u.grps) - 1
    H.s <- lapply(2:length(u.grps), function(j) {
        Xj <- rhs[, grps %in% u.grps[1:j]]
        qrX <- qr(Xj, tol = TOL)
        Q <- qr.Q(qrX)
        tcrossprod(Q[, 1:qrX$rank])
    })
    if (inherits(lhs, "dist")) {
        if (any(lhs < -TOL)) 
            stop("dissimilarities must be non-negative")
        dmat <- as.matrix(lhs^2)
    }else {
        dist.lhs <- as.matrix(vegdist(lhs, method = method))
        dmat <- dist.lhs^2
	   write.table(dist.lhs,"sample.dist.table",sep="\t",quote=F) #########output the distance matrix#############
    }

    n <- nrow(dmat)
    I <- diag(n)
    ones <- matrix(1, nrow = n)
    A <- -(dmat)/2
    G <- -0.5 * dmat %*% (I - ones %*% t(ones)/n)
    SS.Exp.comb <- sapply(H.s, function(hat) sum(G * t(hat)))
    SS.Exp.each <- c(SS.Exp.comb - c(0, SS.Exp.comb[-nterms]))
    H.snterm <- H.s[[nterms]]
    if (length(H.s) > 1) 
        for (i in length(H.s):2) H.s[[i]] <- H.s[[i]] - H.s[[i - 
            1]]
    SS.Res <- sum(G * t(I - H.snterm))
    df.Exp <- sapply(u.grps[-1], function(i) sum(grps == i))
    df.Res <- n - qrhs$rank
    if (inherits(lhs, "dist")) {
        beta.sites <- qr.coef(qrhs, as.matrix(lhs))
        beta.spp <- NULL
    }else {
        beta.sites <- qr.coef(qrhs, dist.lhs)
        beta.spp <- qr.coef(qrhs, as.matrix(lhs))
    }
    colnames(beta.spp) <- colnames(lhs)
    colnames(beta.sites) <- rownames(lhs)
    F.Mod <- (SS.Exp.each/df.Exp)/(SS.Res/df.Res)
    f.test <- function(tH, G, df.Exp, df.Res, tIH.snterm) {
        (sum(G * tH)/df.Exp)/(sum(G * tIH.snterm)/df.Res)
    }
    SS.perms <- function(H, G, I) {
        c(SS.Exp.p = sum(G * t(H)), S.Res.p = sum(G * t(I - H)))
    }
    if (missing(strata))
        strata <- NULL

permuted.index<-function (n, strata)
{
    if (missing(strata) || is.null(strata))
        out <- sample(n, n)
    else {
        out <- 1:n
        inds <- names(table(strata))
        for (is in inds) {
            gr <- out[strata == is]
            if (length(gr) > 1) 
                out[gr] <- sample(gr, length(gr))
        }
    }
    out
}


    p <- sapply(1:permutations, function(x) permuted.index(n,
        strata = strata))
    tH.s <- sapply(H.s, t)
    tIH.snterm <- t(I - H.snterm)
    f.perms <- sapply(1:nterms, function(i) {
        sapply(1:permutations, function(j) {
            f.test(H.s[[i]], G[p[, j], p[, j]], df.Exp[i], df.Res,
                tIH.snterm)
        })
    })
    f.perms <- round(f.perms, 12)
    F.Mod <- round(F.Mod, 12)
    SumsOfSqs = c(SS.Exp.each, SS.Res, sum(SS.Exp.each) + SS.Res)
    tab <- data.frame(Df = c(df.Exp, df.Res, n - 1), SumsOfSqs = SumsOfSqs,
        MeanSqs = c(SS.Exp.each/df.Exp, SS.Res/df.Res, NA), F.Model = c(F.Mod,
            NA, NA), R2 = SumsOfSqs/SumsOfSqs[length(SumsOfSqs)],
        P = c((rowSums(t(f.perms) >= F.Mod) + 1)/(permutations +
            1), NA, NA))
    rownames(tab) <- c(attr(attr(rhs.frame, "terms"), "term.labels")[u.grps],
        "Residuals", "Total")
    colnames(tab)[ncol(tab)] <- "Pr(>F)"
    class(tab) <- c("anova", class(tab))
	tab
}


#.libPaths("/home/zhangjf/R/x86_64-unknown-linux-gnu-library/2.15")
library(vegan)

X <- read.table(argv[2])
X <- data.frame(X)

sampledist <- read.table(argv[1])
#sampledist <- as.dist(sampledist)
dim(X)
lab <- colnames(X)
prefix <- as.vector(argv[3])

colname = c("phenotype","Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")
write.table(t(colname),paste("perm_",prefix,"_","single.txt",sep=""),quote=F,sep="\t",append=T,col.names=F,row.names=F)
for(i in 1:ncol(X)){
	set.seed(0)
    m=which(!is.na(X[,i]))
    sampledist_tem=as.dist(sampledist[m,m])
    if(length( levels(as.factor(X[m,i])))==1){next}
    run=paste("X_tem=data.frame(",lab[i],"=X[m,i])",sep="")
    eval(parse(text = run))
	run <- paste("adonis1(","sampledist_tem~",lab[i],",data=X_tem,","permutations=9999)",sep="")
	tab = eval(parse(text = run))
	tabb <- tab[1,]
	names(tabb) = NULL
	write.table(tabb,paste("perm_",prefix,"_","single.txt",sep=""),quote=F,sep="\t",append=T)
}

#colname = c("phenotype","Df","SumsOfSqs","MeanSqs","F.Model","R2","Pr(>F)")
#write.table(t(colname),paste("perm_",prefix,"_","adjust.txt",sep=""),quote=F,sep="\t",append=T,col.names=F,row.names=F)
#for(i in 1:ncol(X)){
#	set.seed(0)
#	run1 <- paste(lab[-i],collapse="+")
#	run2 <- paste(run1,"+",lab[i],sep="")
#	run <- paste("adonis1(","sampledist~",run2,",data=X,","permutations=9999)",sep="")
#	tab = eval(parse(text = run))
#	ds <- dim(tab)
#	tabb = tab[(ds[1]-2),]
#	names(tabb) = NULL
#	write.table(tabb,paste("perm_",prefix,"_","adjust.txt",sep=""),quote=F,sep="\t",append=T)
#}


