mapplot <- function(){
    par(mar = rep(0, 4),
        oma = rep(0, 4))

    xlim <- c(0, 10)
    ylim <-c(5.05, 10.25)

    plot(0, type = "n", axes = FALSE,
         xlab = "", ylab = "",
         xlim = xlim,    
         ylim = ylim)

    y <- 10

    g <- c(0, 10)

    exons <- list(c(0, 2),
                  c(4.5, 5),
                  c(6, 7.1),
                  c(9, 10))

    trans <- list(c(1, 2, 3, 4),
                  c(1, 3, 4),
                  c(1, 4))

    prot <- list(exons[[1]] + c(0.5, 0), ## 5'UTR
                 exons[[3]],
                 exons[[4]] + c(0, -0.5)) ## 3'UTR

    peps <- list(c(1, 1.29),
                 c(6.1, 6.3), c(6.85, 7),
                 c(9.1, 9.25))


    rect(g[1], y-0.25, g[2], y+0.25, col = "black")
    abline(h = y, lwd = 2)
    text(g[1], y+0.25, expression(g[s]), pos = 3)
    text(g[2], y+0.25, expression(g[e]), pos = 3)
    text(0-0.26, y, "G", pos = 3)

    abline(v = c(unlist(exons),
               unlist(prot),
               unlist(peps)),
           lty = "dotted",
           col = "lightgrey")

    y <- y - 1
    for (i in 1:length(trans)) {
        tr <- trans[[i]]
        abline(h = y, lty = "dotted")
        ## text(0-0.26, y, expression(T[x]), pos = 3)
        text(0-0.26, y,
             substitute(paste("T", list(x)), list(x = i)),
             pos = 3)
        for (j in 1:length(tr)) {
            e1 <- exons[[tr[j]]][1]
            e2 <- exons[[tr[j]]][2]
            rect(e1, y-0.25, e2, y+0.25, col = "grey")
            if (i == 1) {
                text(e1, y+0.25, expression(e[s]^i), pos = 3)
                text(e2, y+0.25, expression(e[e]^i), pos = 3)
            }
            text(mean(c(e1, e2)), y, paste0("i = ", tr[j]), cex = .7)
        }
        y <- y - .6
    }

    y <- y - .4

    abline(h = y, lty = "dotted")
    text(0-0.26, y, expression(P), pos = 3)
    for (i in 1:length(prot)) {
        pr <- prot[[i]]
        rect(pr[1], y - 0.25,
             pr[2], y + 0.25,
             col = "steelblue")
        text(pr[1], y + 0.25, expression(p[s]^j), pos = 3)
        text(pr[2], y + 0.25, expression(p[e]^j), pos = 3)
        text(mean(c(pr[1], pr[2])), y, paste0("j = ", i), cex = .7)
    }
    y <- y - .75

    ## concat protein
    ## center
    x <- mean(xlim)
    protlen <- sum(sapply(prot, diff))
    X0 <- x0 <- x - protlen/2 ## left start

    abline(h = y, lty = "dotted")
    text(0-0.26, y, expression(P), pos = 3)

    for (i in 1:length(prot)) {
        pr <- prot[[i]]
        .pr1 <- x0
        .pr2 <- x0 + (pr[2] - pr[1])
        rect(.pr1, y - 0.25,
             .pr2, y + 0.25,
             col = "steelblue")
        segments(.pr1, y + 0.25, pr[1], y + .75 - .25, lty = "dotted")
        segments(.pr2, y + 0.25, pr[2], y + .75 - .25, lty = "dotted")
        text(mean(c(.pr1,.pr2)), y, paste0("j = ", i), cex = .7)
        x0 <- .pr2
    }

    text(X0, y, expression(1), pos = 2)
    text(X0 + protlen, y, expression(L[P]), pos = 4)

    ## pos and length
    relpeps <- list(c(0.5, 0.29),
                    ## 1.% is length of exon 1
                    c(1.5 + 0.1, 0.2),
                    c(1.5 + 0.85, 0.15),
                    ## 1.1 is length of exon 2
                    c(1.5 + 1.1 + 0.1, 0.15))

    for (i in 1:length(relpeps)) {
        rp <- relpeps[[i]]
        pep <- peps[[i]]
        rect(X0 + rp[1], y - 0.25,
             X0 + rp[1] + rp[2], y + 0.25,
             col = "#FFA50450", lwd = 0)
        segments(X0 + rp[1], y - 0.25,
                 pep[1], y - .75 + .25, lty = "dotted")
        segments(X0 + rp[1] + rp[2], y - 0.25,
                 pep[2], y - .75 + .25, lty = "dotted")
    }
    y <- y - .75


    abline(h = y, lty = "dotted")
    text(0-0.26, y, expression(Pi), pos = 3)
    for (i in 1:length(peps)) {
        pep <- peps[[i]]
        rect(pep[1], y - 0.25,
             pep[2], y + 0.25,
             col = "#FFA504FF")
        text(pep[1], y + 0.25, expression(pi[s]^k), pos = 3, cex = .7)
        text(pep[2], y + 0.25, expression(pi[e]^k), pos = 3, cex = .7)
        text(mean(c(pep[1], pep[2])), y-0.35, paste0("k = ", i), cex = .8)
    }
}

pplot <- function(){
    par(mar = rep(0, 4),
        oma = rep(0, 4))
    ylim <- xlim <- c(0, 10)
    plot(0, type = "n", axes = FALSE,
         xlab = "", ylab = "",
         xlim = xlim,    
         ylim = ylim)

    for (i in 0:9) {
        rect(2, i, 5, i+0.97, col = "lightgrey")
        text(2, i+0.5, paste0("Protein ", 10-i),
             pos = 4)
    }

    rect(0, 0, 1.8, 9.97, lwd = 4)

    for (i in 0:8)
        segments(0, i+0.985, 1.8, i+0.985,
                 lty = "dotted")

    darrow <- function(x1, y1, x2, y2, ...) {
        points(x1, y1, pch = 19)
        arrows(x1, y1, x2, y2, ...)
    }

    for (i in c(0, 1, 3, 4, 5, 7, 8))
        darrow(3.5, i+.5, 6, i+.5, lty = "dotted")

    darrow(3.5, 9.5, 6.5, 9.5)
    darrow(3.5, 6.5, 6.5, 6.5)
    darrow(3.5, 2.5, 6.5, 2.5)

    addpeps <- function(i, j) {
        if (missing(j))
            j <- i + 2
        rect(6.5, i, 8.5, j, col = "lightgrey")
        for (k in seq((i+.25), (j-.25), 0.25))
            segments(6.5, k, 8.5, k)
        text(6.5, k+0.125, expression(peptide[1]), pos = 4, cex = .7)
        text(6.5, k-0.125, expression(peptide[2]), pos = 4, cex = .7)

        rect(8.7, i, 10, j, lwd = 3)
        for (k in seq((i+.25), (j-0.25), 0.25))
            segments(8.7, k, 10, k, lty = "dotted")
    }

    addpeps(7.75)

    addpeps(4.75)

    addpeps(0.75)
}
