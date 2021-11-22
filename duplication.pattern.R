# Is there phylogenetic pattern to gene duplication and loss?

{ # Load packages
  library(phylosignal)
  library(adephylo)
  library(ape)
  library(phylobase)
}

tree <- read.tree("organismal.relationships.tre")
tree <- compute.brlen(tree, method = 1)
plot(tree, use.edge.length = TRUE)

# Duplications indicated in Figure 2, excluding genes limited to
# particular orders, including sis-a, upd. Also excluding tra and the
# "miscellaneous" genes
{
  dup <- rep(0,length(tree$tip.label))
  names(dup) <- tree$tip.label
  # Excluding taxa with unreliable genomic annotations
  # phyloSignal can't handle NA's - They will be left as zero counts.
  # dup[c("Zygentoma","Acridoidea","Grylloidea","Gerroidea","Reduvioidea","Miroidea","Coreoidea","Platypezoidea")] <- NA
  dup[c("Ephydroidea","Papilionoidea","Staphylinoidea","Cucujoidea",
        "Pentatomoidea","Membracoidea","Aphidoidea")] <- 1
  dup[c("Lygaeoidea","Fulgoroidea")] <- 2
  dup[c("Psylloidea","outgroup")] <- 3

  loss <- rep(0,length(tree$tip.label))
  names(loss) <- tree$tip.label
  # Excluding taxa with unreliable genomic annotations
  # phyloSignal can't handle NA's - They will be left as zero counts.
  # loss[c("Zygentoma","Acridoidea","Grylloidea","Gerroidea","Reduvioidea","Miroidea","Coreoidea","Platypezoidea")] <- NA
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
#                  Cmean          I          K    K.star    Lambda
# duplication 0.09900048 0.007968828 1.1813381 0.3772649 0.9867202
# loss        0.22917464 0.029613784 0.3282214 0.3893760 0.7022245
#
# $pvalue
#             Cmean      I      K  K.star  Lambda
# duplication 0.093  0.063  0.010*  0.042*  0.001*
# loss        0.005* 0.010* 0.109   0.019*  0.001*

# Locating the signal with LIPA
# Computes Local Indicator of Phylogenetic Association (local Moran's I) for each tip
lipa <- lipaMoran(p4d, prox.phylo = "patristic", reps = 1e4-1)
{
  barplot.phylo4d(p4d, tree.ladderize = TRUE, use.edge.length = FALSE,
                  bar.col = ifelse(lipa$p.value < 0.05, autumn.maple[5], autumn.maple[1]),
                  trait.bg.col = "gray95",
                  tip.font = 1, center = FALSE , scale = FALSE)
  # Node labels
  text(x=5.5, y=41.5, labels = "Hemiptera",
       col = autumn.maple[5], font = 2, adj = 1, cex = 0.85)
  text(x=6.5, y=30.75, labels = "Hymenoptera",
       col = autumn.maple[1], font = 2, adj = 1, cex = 0.85)
  text(x=6.25, y=20.5, labels = "Coleoptera",
       col = autumn.maple[5], font = 2, adj = 1, cex = 0.85)
}
# Export images as "FigSX.local.phylo.signal.pdf", 8" x 8"

# Local Moran's I and p-values for significant taxa
# duplication
lipa$lipa[which(lipa$p.value[,1]<0.05),1]
# Pentatomoidea    Lygaeoidea  Membracoidea   Fulgoroidea    Psylloidea    Aphidoidea
#    0.12584738    0.05594779    0.20636769    0.22176961    0.15081580    0.27518605
lipa$p.value[which(lipa$p.value[,1]<0.05),1]
# Pentatomoidea    Lygaeoidea  Membracoidea   Fulgoroidea    Psylloidea    Aphidoidea
#        0.0170        0.0375        0.0089        0.0084        0.0138        0.0056
# loss
lipa$lipa[which(lipa$p.value[,2]<0.05),2]
# Tenebrionoidea     Cucujoidea Curculionoidea Chrysomeloidea Elateroidea   Buprestoidea
#      0.2222085      0.2765208      0.3107584      0.3107584   0.1551682      0.1551682
lipa$p.value[which(lipa$p.value[,2]<0.05),2]
# Tenebrionoidea     Cucujoidea Curculionoidea Chrysomeloidea Elateroidea   Buprestoidea
#         0.0016         0.0012         0.0021         0.0020      0.0102         0.0105

