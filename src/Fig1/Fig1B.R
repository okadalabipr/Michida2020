library(venneuler)
#SE
n.A = 147
n.B = 20
n.AB = 113

v <- venneuler(c(A = n.A, B = n.B , "A&B" = n.AB))
v$labels <- factor(c("", ""))
png("Venn_SE.png", w=1200, h=1200)
plot(v, col = c("deeppink", "deepskyblue2", alpha = 0.5))

dev.off()

n.A = 4812
n.B = 1804
n.AB = 2260

#TE
v <- venneuler(c(A = n.A, B = n.B , "A&B" = n.AB))
v$labels <- factor(c("", ""))
png("Venn_TE.png", w=1200, h=1200)
plot(v, col = c("deeppink", "deepskyblue2", alpha = 0.5))

dev.off()

#Obtained figures were changed their aspect ration in the same ratio by PowerPoint (MicroSoft).