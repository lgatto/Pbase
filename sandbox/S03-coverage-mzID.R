library("mzID")
library("Biostrings")
library("Pbase")

system.time(
p <- Proteins("swissprot_human_canonical_19_09_12.fasta")
)
system.time(
m <- mzID("Thermo_Hela_PRTC_1.mzid")
)
system.time(
  pc <- proteinCoverage(p, m)
)

cv <- acols(pc)$Coverage
n <- elementLengths(aa(pc))

par(mfrow = c(2, 2))
boxplot(cv)
plot(n, cv, col = "#00000040")
boxplot(cv[cv > 0])
plot(n[cv > 0], cv[cv > 0], col = "#00000040")
par(mfrow = c(1, 1))

pcols(pc[1:5])
