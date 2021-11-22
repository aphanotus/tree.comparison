# Is there phylogenetic pattern to gene duplication and loss?

{ # Load packages
  library(phylosignal)
  library(adephylo)
  library(ape)
  library(phylobase)
  library(dplyr)
  library(magrittr)
}

tree <- read.tree("organismal.relationships.tre")

# Excluding taxa with unreliable genomic annotations
tree <- drop.tip(tree,c("Zygentoma","Acridoidea","Grylloidea","Gerroidea","Reduvioidea","Miroidea","Coreoidea","Platypezoidea"))
tree <- compute.brlen(tree, method = 1)
plot(tree, use.edge.length = TRUE)

# Duplications indicated in Figure 2, excluding genes limited to
# particular orders, including sis-a, upd. Also excluding tra and the
# "miscellaneous" genes
{
  dup <- rep(0,length(tree$tip.label))
  names(dup) <- tree$tip.label
  dup[c("Ephydroidea","Papilionoidea","Staphylinoidea","Cucujoidea",
        "Pentatomoidea","Membracoidea","Aphidoidea")] <- 1
  dup[c("Lygaeoidea","Fulgoroidea")] <- 2
  dup[c("Psylloidea","outgroup")] <- 3

  loss <- rep(0,length(tree$tip.label))
  names(loss) <- tree$tip.label
  loss[c("Thysanoptera","Aleyrodoidea","Aphidoidea",
         "Cimicoidea","Siphonaptera","Muscoidea")] <- 1
  loss[c("Tenthredinoidea","Orussoidea","Chalcidoidea","Vespoidea","Apoidea","Formicoidea",
         "Yponomeutoidea")] <- 1
  loss[c("Ephemeroptera","Odonata","Cimicoidea","Pentatomoidea","Lygaeoidea","Phthiraptera",
         "Cephoidea")] <- 2
  loss[c("Buprestoidea","Elateroidea","Scarabaeoidea","Tenebrionoidea","Cucujoidea",
         "Chrysomeloidea","Curculionoidea")] <- 2
  loss[c("Membracoidea")] <- 3

  dat <- data.frame(duplication = dup, loss = loss)
}

# See http://www.francoiskeck.fr/phylosignal/demo_general.html
p4d <- phylo4d(tree, dat)
phyloSignal(p4d = p4d, method = "all")
# $stat
#                 Cmean          I         K    K.star   Lambda
# duplication 0.1512521 0.02647819 1.2288359 0.4376425 1.033838
# loss        0.1512521 0.02647819 1.2288359 0.4376425 1.033838
#
# $pvalue
#             Cmean      I      K  K.star  Lambda
# duplication 0.072  0.017* 0.011*  0.021*  0.001*
# loss        0.001* 0.002* 0.039*  0.001*  0.001*

# Locating the signal with LIPA
# Computes Local Indicator of Phylogenetic Association (local Moran's I) for each tip
lipa <- lipaMoran(p4d, prox.phylo = "patristic", reps = 1e4-1)

dat$dup.I <- NA; dat$dup.p <- NA; dat$loss.I <- NA; dat$loss.p <- NA
for (i in 1:(dim(dat)[1])) {
  taxon.i <- rownames(dat)[i]
  dat$dup.I[i] <- lipa$lipa[taxon.i,1]
  dat$dup.p[i] <- lipa$p.value[taxon.i,1]
  dat$loss.I[i] <- lipa$lipa[taxon.i,2]
  dat$loss.p[i] <- lipa$p.value[taxon.i,2]
}

# Local Moran's I and p-values for significant taxa
# duplication
filter(dat, duplication > 0, dup.p < 0.05) %>% select(dup.I, dup.p)
#                   dup.I  dup.p
# Lygaeoidea    0.2998664 0.0054
# Pentatomoidea 0.2342867 0.0085
# Membracoidea  0.2334781 0.0095
# Fulgoroidea   0.3151001 0.0052
# Psylloidea    0.2233402 0.0082
# Aphidoidea    0.2756377 0.0050

# loss
filter(dat, loss > 0, loss.p < 0.05) %>% select(loss.I, loss.p)
#                    loss.I loss.p
# Curculionoidea 0.24214131 0.0070
# Chrysomeloidea 0.24214131 0.0073
# Cucujoidea     0.20819483 0.0047
# Tenebrionoidea 0.15534267 0.0076
# Elateroidea    0.09463598 0.0455
# Buprestoidea   0.09463598 0.0459
# Lygaeoidea     0.21701026 0.0102
# Pentatomoidea  0.21701026 0.0101
# Cimicoidea     0.16957149 0.0129

{
  barplot.phylo4d(p4d, tree.ladderize = TRUE, use.edge.length = FALSE,
                  bar.col = ifelse(lipa$p.value < 0.05, autumn.maple[5], autumn.maple[1]),
                  trait.bg.col = "gray95",
                  tip.font = 1, center = FALSE , scale = FALSE)
  # Node labels
  text(x=4.5, y=36.75, labels = "Hemiptera",
       col = autumn.maple[5], font = 2, adj = 1, cex = 0.85)
  text(x=4.75, y=29.6, labels = "Hymenoptera",
       col = autumn.maple[1], font = 2, adj = 1, cex = 0.85)
  text(x=4.5, y=19.5, labels = "Coleoptera",
       col = autumn.maple[5], font = 2, adj = 1, cex = 0.85)
}
# Export images as "FigS2.local.phylo.signal.pdf", 8" x 8"
