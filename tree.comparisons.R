# Compare the superfamily-level topology of each gene tree to a
# consensus organismal tree (Misof et al. 2014) to find the
# tree distance as normalized Clustering Information Distance
# (Smith 2020 https://doi.org/10.1093/bioinformatics/btaa614 )

{ # Load packages
  library(tidyverse)
  library(phytools)
  library(phylotools)
  library(ips)
  library(TreeDist)
  library(paletteer)
  library(ggpubr)
  library(ggrepel)
  library(ggtree)
}

# All trees must be in Newick format -- not Nexus
# These are expected to be "RAxML_bipartitions." files

##########################
### Functions
##########################

quantify.polytomy <- function(tree) {
  tree <- unroot(tree)
  return(1 - ( (tree$Nnode -1) / (length(tree$tip.label)-3) ))
}

remove.redundant.parentheses <- function(s) {
  # Remove parentheses wrapping a single taxon
  continue <- grepl("\\(\\w+\\)",s)
  while (continue) {
    s <- gsub("\\((\\w+)\\)", "\\1", s)
    continue <- grepl("\\(\\w+\\)",s)
  }
  sv <- unlist(lapply(1:nchar(s), function(x) substr(s,x,x)))
  if (length(grep("\\(",sv))<2) {
    return(paste0(sv, collapse = ""))
  } else {
    parentheses.matrix <- data.frame(matrix(nrow=length(grep("\\(",sv)), ncol=2,
                                            dimnames=list(1:length(grep("\\(",sv)),c("open","close"))))
    matrix.pos <- 1
    for (i in 1:(length(sv)-1)) {
      if (sv[i]=="(") {
        parentheses.matrix[matrix.pos,1] <- i
        counter <- 1
        j <- i
        while (counter > 0) {
          j <- j+1
          if (sv[j]=="(") { counter <- counter +1 } else {
            if (sv[j]==")") { counter <- counter -1 }
          }
          parentheses.matrix[matrix.pos,2] <- j
        }
        matrix.pos <- matrix.pos +1
      }
    }
    parantheses.to.cut <- NULL
    for (i in 1:(dim(parentheses.matrix)[1]-1)) {
      if (parentheses.matrix$open[i] == parentheses.matrix$open[i+1]-1) {
        if (parentheses.matrix$close[i] == parentheses.matrix$close[i+1]+1)
          parantheses.to.cut <- c(parantheses.to.cut, parentheses.matrix$open[i], parentheses.matrix$close[i])
      }
    }
    if (is.null(parantheses.to.cut)) { return(paste0(sv, collapse = "")) }
    else { return(paste0(sv[-parantheses.to.cut], collapse = "")) }
  }
} # End function  remove.redundant.parentheses

vectorize.newick.tips <- function(file = "tmp.tre",
                                  sampled.taxonomy,
                                  include.topology = TRUE,
                                  output.display.taxa = TRUE) {
  # Read in the tree file as a string in order to substitute taxon names
  s <- read.delim(file, sep =";", header = FALSE)[1]
  # Remove any branch lengths
  s <- gsub("\\:\\d+\\.\\d+","",s)
  # Remove node supports
  s <- gsub("\\)\\d+",")",s)
  # Split the string into a vector
  s <- unlist(str_split(s, ","))
  # # Change Drosophila abbreviation to full genus name
  # s <- sub("D_","Drosophila_",s)
  # s <- sub("Dmel_","Drosophila_",s)
  for (j in 1:length(s)) {
    # Note the number of opening and closing parentheses
    open.parantheses <- str_count(s[j],"\\(")
    close.parantheses <- str_count(s[j],"\\)")
    # Find the matching entry in the taxonomy dataset
    x <- which(unlist(lapply(sampled.taxonomy$genus, function (x) { grepl(pattern=x, x=s[j]) })))[1]
    # if ((length(x) != 1) | is.na(x)) {
    #   cat(paste0(gsub("[\\(\\)]","",s[j])," not found in taxonomy\n"))
    # }
    if (output.display.taxa) {
      if (include.topology) {
        s[j] <- paste0(
          paste0(rep("(",open.parantheses),collapse = ""),
          sampled.taxonomy$display.taxon[x],
          paste0(rep(")",close.parantheses),collapse = "")
        )
      } else {
        s[j] <- sampled.taxonomy$display.taxon[x]
      }
    } else {
      x <- paste0(str_split_fixed(s[j],"_",3)[,1:2], collapse = "_")
      if (include.topology) {
        s[j] <- paste0(
          paste0(rep("(",open.parantheses),collapse = ""),
          x,
          paste0(rep(")",close.parantheses),collapse = "")
        )
      } else {
        s[j] <- gsub("\\(","",x)
      }
    }
  } # End j loop
  return(s)
} # End of function  vectorize.newick.tips

clean.up.newick.vector <- function(s) {
  # Collapse redundancies
  for (j in 1:(length(s)-1)) {
    s.j <- gsub("[()]","",s[j])
    s.j1 <- gsub("[()]","",s[j+1])
    if ( (s.j == s.j1) ) {
      s[j] <- sub(s.j,"",s[j])
    }
  }

  # Reassemble the string
  s <- paste0(s, collapse = ",")
  # Clean up by removing commas adjacent to parentheses
  s <- gsub("\\(\\,","(",s)
  s <- gsub("\\,\\)",")",s)
  continue <- grepl("\\(\\)",s)
  while (continue) {
    s <- gsub("\\(\\)","",s)
    s <- gsub("\\(\\,","(",s)
    s <- gsub("\\,\\)",")",s)
    continue <- (grepl("\\(\\)",s) | grepl("\\(\\,",s) | grepl("\\,\\)",s))
  }

  # Remove parentheses wrapping a single taxon
  continue <- grepl("\\(\\w+\\)",s)
  while (continue) {
    s <- gsub("\\((\\w+)\\)", "\\1", s)
    continue <- grepl("\\(\\w+\\)",s)
  }

  # Remove any double commas
  continue <- grepl("\\,\\,",s)
  while (continue) {
    s <- gsub("\\,\\,", ",", s)
    continue <- grepl("\\,\\,",s)
  }

  # Remove commas in invalid positions
  s <- gsub("\\(\\,", "(", s)
  s <- gsub("\\,\\)", ")", s)

  return(s)
} # End of function  clean.up.newick.vector

compare.gene.trees.to.organismal.tree <- function (
  input.tree.folder = "./analysis.trees/",
  organismal.tree.file = "organismal.relationships.tre",
  taxonomy.file = "taxids.csv",
  support.cutoff = 50,
  bootstrap.reps = 1e4,
  condensed.trees.output.folder = "./condensed.trees"
)
{
  # Input files
  tree.files <- paste0(input.tree.folder,list.files(input.tree.folder))
  sampled.taxonomy <- read.csv(taxonomy.file)
  sampled.taxonomy$genus <- str_split_fixed(sampled.taxonomy$species,"[\\s_]",2)[,1]
  organismal.tree <- unroot(read.tree(organismal.tree.file))

  # Get gene names from the tree files
  x <- sub(input.tree.folder,"",tree.files)
  x <- sub("RAxML_bipartitions.","",x)
  gene.names <- str_split_fixed(x,"\\.",2)[,1]

  # Initialize a data frame to capture the results
  x <- c("gene","process","group","species","taxa","polytomy","paraphyly",
         "MCI","MCI.ci", "nCID","nCID.ci" )
  tree.comparisons <- data.frame(matrix(ncol = length(x),
                                        nrow = length(tree.files)))
  colnames(tree.comparisons) <- x

  # Initialize a list to capture distributions from nCID sampling
  nCID.distributions <- list(
    # values = list(),
    plots = list()
  )

  # Main loop to process each tree file
  options(warn=-1)
  par(mfrow=c(1,2))
  for (i in 1:length(tree.files)) {

    tree.comparisons$gene[i] <- gene.names[i]
    tree.i <- read.tree(tree.files[i])

    # Plot the input tree with node supports
    tmp.tree <- tree.i
    tmp.tree$tip.label <- str_split_fixed(tmp.tree$tip.labe,"_",2)[,1]
    plot(tmp.tree, cex = 0.8, main = paste0(gene.names[i]), no.margin = FALSE, use.edge.length = FALSE)
    nodelabels(text = tmp.tree$node.label, cex = 0.7,
               node = 1:tmp.tree$Nnode+Ntip(tmp.tree),
               frame="none", adj=c(1.1,-0.4), col="darkred", font =2)

    # Customize the organismal tree to match taxa in the gene tree
    org.tree.i <- organismal.tree

    # Find analysis level taxa for each tip, and remove
    # any taxa that aren't covered in the organismal tree
    s <- vectorize.newick.tips(file = tree.files[i],
                               sampled.taxonomy = sampled.taxonomy,
                               include.topology = FALSE,
                               output.display.taxa = TRUE)
    x <- which(s %in% org.tree.i$tip.label)
    tree.i <- keep.tip(tree.i, x)

    # Record the number of species / tips in the tree
    tree.comparisons$species[i] <- length(tree.i$tip.label)

    # Collapse poorly supported nodes into polytomies
    tree.i <- ips::collapseUnsupportedEdges(tree.i, cutoff = support.cutoff)

    # Write a temporary tree file, then read in the tree file
    # as a string in order to substitute taxon names
    write.tree(tree.i,"tmp.tre")
    s <- vectorize.newick.tips(file = "tmp.tre",
                               sampled.taxonomy = sampled.taxonomy,
                               include.topology = TRUE,
                               output.display.taxa = TRUE)
    s <- clean.up.newick.vector(s)

    # Check that the tree has more than 1 tip
    if (str_count(s,"\\,") < 2) {
      tree.comparisons$taxa[i] <- 1
    } else {

      # Convert to a phylo (tree) data structure
      s <- remove.redundant.parentheses(s)
      tree.i <- unroot(read.tree(text = paste0(s,";")))

      # Ensure the gene tree only has the taxa in the organismal tree
      x <- which(tree.i$tip.label %in% org.tree.i$tip.label)
      tree.i <- keep.tip(tree.i, x)

      # Record the number of unique taxa at the analytical level (includes outgroup)
      tree.comparisons$taxa[i] <- length(unique(tree.i$tip.label))

      # Remove redundant sister taxa
      {
        # This occurs when multiple species represent one analysis taxon and it's not paraphyletic
        tips.to.remove <- NULL
        for (j in 1:length(tree.i$tip.label)) {
          sister.numbers <- getSisters(tree.i, j, mode = "number")
          sister.numbers <- sister.numbers[which(sister.numbers <= length(tree.i$tip.label) )]
          names(sister.numbers) <- getSisters(tree.i, j, mode = "label")$tips
          if (!(j %in% tips.to.remove)) {
            tips.to.remove <- c(tips.to.remove, sister.numbers[which(names(sister.numbers)==tree.i$tip.label[j])])
          }
        }
        tree.i <- drop.tip(tree.i, tips.to.remove)
      }

      # Ensure the organismal tree only has the taxa in the gene tree
      x <- which(org.tree.i$tip.label %in% tree.i$tip.label)
      org.tree.i <- keep.tip(org.tree.i, x)

      # Save the condensed tree, including all paraphyly
      if (!dir.exists(condensed.trees.output.folder)) { dir.create(condensed.trees.output.folder) }
      write.tree(tree.i, file = paste0(condensed.trees.output.folder,"/",gene.names[i],".topology.bs",support.cutoff,".tre"))

      # Record polytomy
      tree.comparisons$polytomy[i] <- quantify.polytomy(tree.i)

      # Record paraphyly
      duplicated.tips <- duplicated(tree.i$tip.label)
      tree.comparisons$paraphyly[i] <- sum(duplicated.tips) / (0.5*length(unique(tree.i$tip.label)))

      plot(tree.i, cex = 0.8, no.margin = FALSE,
           main = paste0(
             "polytomy: ",signif(tree.comparisons$polytomy[i],2),
             "\nparaphyly: ",signif(tree.comparisons$paraphyly[i],2)) )

      # Check that the tree has at least 4 tips
      if (length(unique(tree.i$tip.label)) < 3) {
        cat(paste0(gene.names[i],":\t",s," only\n"))
      } else { # Make meaningful comparisons of the trees

        # The tree distance calculation requires trees to have unique (monophyletic) tips.
        # So in cases of paraphyly, repeated tip labels are removed, randomly keeping one.
        # Therefore, MCI and nCID are taken as the means calculated from multiple iterations.
        sample.tree.dist <- FALSE
        if (sum(duplicated.tips) > 0) { sample.tree.dist <- TRUE }
        MCI.samples <- vector()
        nCID.samples <- vector()
        # The permutations are also used to find the 95% CI bound using bootstrapping
        MCI.ci <- vector()
        nCID.ci <- vector()
        unique.tip.labels <- unique(tree.i$tip.label)
        tree.i.backup <- tree.i
        org.tree.i.backup <- org.tree.i
        # Permutation loop
        progress <- txtProgressBar(min = 1, max = bootstrap.reps, initial = 1, char = ".", style = 3)
        for (j in 1:bootstrap.reps) {
          setTxtProgressBar(progress, j)

          # If there's no paraphyly, then tree distances are calculated once,
          # after the permutation loop
          if (sample.tree.dist) {
            # Remove paraphyly by randomly retaining one of each repeated taxon
            for (k in 1:length(unique.tip.labels)) {
              tips.k <- which(tree.i$tip.label == unique.tip.labels[k])
              if (length(tips.k) > 1) {
                tips.k <- sample(x = tips.k, size = length(tips.k)-1, replace = FALSE)
                tree.i <- drop.tip(tree.i, tips.k)
              }
            } # End k loop to remove paraphyly

            # Tree distances
            MCI.samples <- c(MCI.samples,MutualClusteringInfo(org.tree.i, tree.i, normalize = FALSE))
            nCID.samples <- c(nCID.samples,ClusteringInfoDistance(org.tree.i, tree.i, normalize = TRUE))
          } # End if(sample.tree.dist)

          # Bootstrapping
          tree.i$tip.label <- sample(tree.i$tip.label, replace = FALSE)
          MCI.ci <- c(MCI.ci,MutualClusteringInfo(org.tree.i, tree.i, normalize = FALSE))
          nCID.ci <- c(nCID.ci,ClusteringInfoDistance(org.tree.i, tree.i, normalize = TRUE))

          # Re-set the trees for the next iteration
          tree.i <- tree.i.backup
          org.tree.i <- org.tree.i.backup

        } # End j loop for permutations
        close(progress)

        # tree distances
        if (sample.tree.dist) {
          tree.comparisons$MCI[i] <- mean(MCI.samples)
          tree.comparisons$nCID[i] <- mean(nCID.samples)

          # Save the distribution of nCID samples
          p <- data.frame(value=nCID.samples) %>%
            ggplot(aes(x=value)) +
            theme_bw() +
            geom_histogram() +
            labs(title = paste0(gene.names[i])) +
            scale_y_continuous(name = NULL) +
            scale_x_continuous(name = "nCID", limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1"))
          x <- length(nCID.distributions$plots)
          # nCID.distributions$values[[x+1]] <- nCID.samples
          # names(nCID.distributions$values[[x+1]]) <- gene.names[i]
          nCID.distributions$plots[[x+1]] <- p
          # names(nCID.distributions$plots[[x+1]]) <- gene.names[i]

        } else {
          tree.comparisons$MCI[i] <- MutualClusteringInfo(org.tree.i, tree.i, normalize = FALSE)
          tree.comparisons$nCID[i] <- ClusteringInfoDistance(org.tree.i, tree.i, normalize = TRUE)
        }

        # Find the 95% (one-sided) quantile of the bootstrap tree distance values
        tree.comparisons$MCI.ci[i] <- quantile(MCI.ci, 0.95, na.rm = TRUE)
        tree.comparisons$nCID.ci[i] <- quantile(nCID.ci, 0.05, na.rm = TRUE)

      } # End else from if (length(unique(tree.i$tip.label)) < 3)

    } # End else from if (str_count(s,"\\,") < 2)

    x <- select(tree.comparisons, -process, -group)
    print(x[i,])

  } # i loop for each gene tree file
  options(warn=0)
  par(mfrow=c(1,1))

  # Delete the tmp tree file
  system2(command = "rm", args = "tmp.tre")

  # Annotate genes based on process
  {
    x <- c("AATS","AcCoA","actb","ADA","ArgK","Arr2","CAD","Cam","Cda4","Cp1","CTR9","CycC","CycE","DDC","DNApol-d","EcR","EF1a1","EF2","eIF3a","eRF1","Fen1","Gel","HDAC1","Iap2","IDH","kkv","mp20","Npl4","PABP1","PEPCK","PGD","PPO1","Prp6","pug","RCC1","RpII18","RpII215","Sf1","SOD1","Spec-a","Syx1A","TPI","UbA5","UDE","Wnt1")
    tree.comparisons$process[which(tree.comparisons$gene %in% x)] <- "reference"
    x <- c("snf","spf45","sisa","scute","upd","runt","dpn","dsx","dsx_OD1","dsx_OD2","ix","fru","Sxl","tra","tra_SDP","tra2","tra2_RRM","vir","fl2d","her", "Masc","chinmo","Mdmd","CWC22","PSI")
    tree.comparisons$process[which(tree.comparisons$gene %in% x)] <- "somatic sex"
    x <- c("otu","ovo","nclb","stil")
    tree.comparisons$process[which(tree.comparisons$gene %in% x)] <- "germline sex"
    x <- c("mle","mof","msl1","msl1_PEHE","msl2","msl2_ZFR","msl3")
    tree.comparisons$process[which(tree.comparisons$gene %in% x)] <- "dosage comp"
    x <- c("dmrt11E","dmrt93B","dmrt99B")
    tree.comparisons$process[which(tree.comparisons$gene %in% x)] <- "other DMRT"
    tree.comparisons <- tree.comparisons %>%
      mutate(process=fct_relevel(tree.comparisons$process, c("reference","other DMRT","dosage comp","germline sex","somatic sex")))
  }

  # Annotate genes into "group" for two-sample stats later.
  # Genes that are secondary to the focus of the project are
  # marked NA to be excluded from those tests.
  tree.comparisons$group <- ifelse((tree.comparisons$process == "reference"),"reference","sex determination")
  tree.comparisons$group[which(grepl("dosage",tree.comparisons$process))] <- NA
  tree.comparisons$group[which(grepl("other",tree.comparisons$process))] <- NA
  tree.comparisons$group[which(grepl("_",tree.comparisons$gene))] <- NA

  # Save an image of nCID sample distributions
  p <- ggarrange(plotlist = nCID.distributions$plots)
  ggsave(paste0("FigSX.nCID.permutations.bs",support.cutoff,".pdf"), p)

  return(tree.comparisons)
} # End of function  compare.gene.trees.to.organismal.tree

# functions to calculate alignment information
parsimony.informative.sites <- function (alignment, ignore.gaps = TRUE, normalize = TRUE) {
  # Parsimony-informative site: contains at least two types of amino acids,
  # and at least two of them occur with a minimum frequency of two.
  # The function takes an alignment in the form provided by phylotools::read.phylip
  # If ignore.gaps = TRUE and normalize = TRUE, eliminates apomorphic sites from the normalization
  if (length(unique(nchar(alignment$seq.text)))!=1) {
    stop("`parsimony.informative.sites` requires aligned sequences.")
  }
  taxa <- length(alignment$seq.text)
  pos <- nchar(alignment$seq.text)[1]
  pos.minus.apomorphy <- pos
  informative.sites <- 0
  for (i in 1:pos) {
    char.i <- substr(alignment$seq.text,i,i)
    if (ignore.gaps) {
      char.i <- char.i[-which(char.i=="-")]
      if (length(char.i)==1) { pos.minus.apomorphy <- pos.minus.apomorphy -1 }
    }
    unique.chars <- length(unique(char.i))
    if (unique.chars>1 & unique.chars<(pos-1)) {
      informative.sites <- informative.sites +1
    }
  }
  if (normalize){
    if (ignore.gaps) {
      return(informative.sites/pos.minus.apomorphy)
    } else {
      return(informative.sites/pos)
    }
  } else {
    return(informative.sites)
  }
} # End function parsimony.informative.sites

mean.pairwise.identity <- function (alignment) {
  # The function takes an alignment in the form provided by phylotools::read.phylip
  # Gaps are not counted as matches
  if (length(unique(nchar(alignment$seq.text)))!=1) {
    stop("`mean.pairwise.identity` requires aligned sequences.")
  }
  taxa <- length(alignment$seq.text)
  pos <- nchar(alignment$seq.text)[1]
  dat <- t(sapply(strsplit(alignment$seq.text,""), as.character))
  id.matrix <- matrix(nrow = taxa, ncol = taxa)
  for (i in 1:(taxa-1)) {
    for (j in (i+1):taxa) {
      numerator <- sum(apply(dat[c(i,j),],2,function(x){ (!all(x=="-")) & (x[1]==x[2])}))
      denominator <- sum(apply(dat[c(i,j),],2,function(x){any(x!="-")}))
      id.matrix[i,j] <- numerator / denominator
    }
  }
  return(mean(id.matrix, na.rm = TRUE))
} # End function mean.pairwise.identity

##########################
### Alignment metrics
##########################
{
  phy.files <- list.files("./analysis.phy/")
  phy.files <- phy.files[grepl("\\.phy$",phy.files)]

  # Get gene names from the phylip files
  phy.names <- sub("blastp\\.\\w+\\.","",phy.files)
  phy.names <- str_split_fixed(phy.names,"\\.",2)[,1]

  seq.identity <- vector()
  info.sites <- vector()
  for (i in 1:length(phy.files)) {
    x <- read.phylip(paste0("analysis.phy/",phy.files[i]))
    seq.identity[i] <- mean.pairwise.identity(x)
    info.sites[i] <- parsimony.informative.sites(x)
    cat(i,phy.names[i],"\tseq identity:",seq.identity[i],"\tinformative sites:",info.sites[i],"\n")
  }
  names(seq.identity) <- phy.names
  names(info.sites) <- phy.names
}

##########################
### Bootstrap cutoff of 50
##########################
# Focus analysis on a bootstrap cutoff of 50, but run at cutoffs of 20 and 90
# for supplemental figures and stats
tree.comp50 <- compare.gene.trees.to.organismal.tree(support.cutoff = 50)

# Combine results
tree.comp50$seq.identity <- NA
tree.comp50$info.sites <- NA
for (i in 1:(dim(tree.comp50)[1])) {
  gene.i <- tree.comp50$gene[i]
  tree.comp50$seq.identity[i] <- seq.identity[gene.i]
  tree.comp50$info.sites[i] <- info.sites[gene.i]
}

# Write the results to a file
(s <- paste0("tree.comparisons.bs50.",format(Sys.time(), "%Y.%m.%d"),".csv"))
write.csv(tree.comp50, s, row.names = FALSE)

# file.choose()
# tree.comp50 <- read.csv("tree.comparisons.bs50.2021.11.11.csv")
# tree.comp50 <- tree.comp50 %>% mutate(process=fct_relevel(tree.comp50$process, "reference"))

# Prepare figures
autumn.maple <-
  c(paletteer_d("dichromat::BluetoGray_8")[7],
    paletteer_d("jcolors::pal12")[c(8,10,12)],
    paletteer_d("dichromat::DarkRedtoBlue_12")[12])
# scales::show_col(autumn.maple)
# scales::show_col(dichromat(autumn.maple))

# Figure with basic stats from each gene tree
df <- tree.comp50 %>%
  select(gene, process, species, taxa, polytomy, paraphyly, seq.identity, info.sites, nCID, nCID.ci ) %>%
  mutate(nCID.by.ID = nCID/seq.identity ) %>%
  arrange(process, desc(nCID.by.ID)) %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  mutate(gene = factor(gene, levels=gene)) %>%
  mutate(sample.size = paste0("n = ",taxa," (",species," sp.)") ) %>%
  pivot_longer(cols = c(5:9,11))

df$nCID.ci[which(df$name != "nCID")] <- NA
df$sample.size[which(df$name != "seq.identity")] <- NA

df <- df %>%
  mutate(name = recode(as.character(name), seq.identity = "sequence identity")) %>%
  mutate(name = recode(as.character(name), info.sites = "informative sites")) %>%
  mutate(name = recode(as.character(name), nCID = "tree distance (nCID)")) %>%
  mutate(name = recode(as.character(name), nCID.by.ID = "nCID / identity")) %>%
  mutate(name = factor(name, levels=unique(name)[c(3,4,1,2,5,6)]))

gene.tree.summary.stats.plot <- df %>%
  ggplot( aes(x = gene, y = value, fill = process) ) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.ticks = element_blank(),
        strip.text = element_text(face="bold"),
        strip.background = element_blank()) +
  facet_grid(.~name, scales = "free_x") +
  geom_bar(stat="identity") +
  scale_fill_manual(values = autumn.maple) +
  guides(fill = guide_legend(reverse=TRUE)) +
  geom_errorbar(aes(ymin=nCID.ci, ymax=nCID.ci),
                color = "grey40", size = 1, width = 0.7) +
  geom_text(aes(x = gene, y = 1, label = sample.size),
            color = "grey65", size = 3, hjust = 1) +
  labs(x = NULL, y = NULL) +
  coord_flip()

gene.tree.summary.stats.plot

ggsave("FigS4.gene.tree.summary.stats.plot.pdf",gene.tree.summary.stats.plot,
       width = 6.5, height = 9, scale = 1.5)

genes.to.label <-
  c("snf","spf45","scute","runt","dpn","dsx","dsx (OD1)","dsx (OD2)","ix","fru","Sxl","tra","tra (SDP)","tra2","vir","fl2d",
    "otu","ovo","nclb","chinmo",
    "mle","mof","msl1","msl2","msl3",
    "dmrt11E","dmrt93B","dmrt99B",
    "Iap2","mp20","RpII18","RpII215","Cam","EF1a1","Wnt1","RCC1","Gel","CAD","Spec-a")

nCIDbyID.vs.info.plot <- tree.comp50 %>%
  mutate(nCID.by.ID = nCID/seq.identity ) %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  mutate(name = ifelse(gene %in% genes.to.label,gene,NA) ) %>%
  ggplot( aes(x = info.sites, y = nCID.by.ID, color = process, label = name) ) +
  theme_bw() +
  theme(legend.justification=c(0,1), legend.position=c(0.01,0.99),
        axis.ticks = element_blank() ) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(size=3, color="grey15", max.iter = 1e6, box.padding = 0.45, max.overlaps = 100) +
  scale_color_manual(values = autumn.maple) +
  guides(color = guide_legend(reverse=TRUE)) +
  ylim(c(0,4)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  labs(x = "portion of parsimony informative sites", y = "tree distance / portion of identical sites")

nCIDbyID.boxplot <- tree.comp50 %>%
  filter(!grepl("_",gene)) %>%
  mutate(nCID.by.ID = nCID/seq.identity ) %>%
  mutate(name = ifelse(gene %in% genes.to.label,gene,NA) ) %>%
  ggplot( aes(x = group, y = nCID.by.ID, color = process, label = name) ) +
  theme_bw() +
  theme(legend.position="none",
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank() ) +
  geom_boxplot(color = "gray50", alpha = 0.85) +
  geom_jitter(height = 0, width = 0.2, size = 3, alpha = 0.85) +
  geom_text(size = 3, color="grey15", hjust = 0, nudge_x = -0.35) +
  scale_color_manual(values = autumn.maple) +
  scale_y_continuous(name = NULL, limits = c(0,4)) + # "tree distance / portion of identical sites"
  scale_x_discrete(name = " ", labels = c("reference","sex determination"," "))

gene.comp.plot <- ggarrange(nCIDbyID.vs.info.plot, nCIDbyID.boxplot, widths=c(2,1))
gene.comp.plot

ggsave("Fig2.genetree.comp.plots.pdf",gene.comp.plot,
       width = 9, height = 6, scale = 1)

# Other supplementary plots:
# - Scatterplot: info.sites vs seq.identity
# - Scatterplot: polytomy vs. vs info.site
# - Scatterplot: nCID vs polytomy
# - Scatterplot: nCID vs info.site

correlation.report <- function(x,y) {
  Rsq <- cor(x,y)^2
  pvalue <- cor.test(x,y)$p.value
  paste0("R^2 = ",signif(Rsq,2),"\np = ",signif(pvalue,2))
}

infoSites.vs.identity.plot.txt <- with(tree.comp50, correlation.report(seq.identity, info.sites))
infoSites.vs.identity.plot <- tree.comp50 %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  mutate(name = ifelse(gene %in% genes.to.label,gene,NA) ) %>%
  ggplot( aes(x = seq.identity, y = info.sites, color = process, label = name) ) +
  theme_bw() +
  theme(#legend.justification=c(0,1), legend.position=c(0.01,0.99),
        axis.ticks = element_blank() ) +
  geom_smooth(method=lm, se=FALSE, color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  annotate(geom="text", x=0.25/2, y= 0.25/2, label = infoSites.vs.identity.plot.txt) +
  geom_text_repel(size=2, color="grey15", max.iter = 1e6, box.padding = 0.45, max.overlaps = 100) +
  scale_color_manual(values = autumn.maple) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  coord_fixed() +
  labs(x = "portion of identical sites", y = "\nportion of parsimony informative sites")

polytomy.vs.info.plot.txt <- with(tree.comp50, correlation.report(info.sites, polytomy))
polytomy.vs.info.plot <- tree.comp50 %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  mutate(name = ifelse(gene %in% genes.to.label,gene,NA) ) %>%
  ggplot( aes(x = info.sites, y = polytomy, color = process, label = name) ) +
  theme_bw() +
  theme(#legend.justification=c(0,0), legend.position=c(0.01,0.01),
        axis.ticks = element_blank() ) +
  # geom_smooth(method=lm, se=FALSE, color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  annotate(geom="text", x=0.25/2, y= 0.25/2, label = polytomy.vs.info.plot.txt) +
  geom_text_repel(size=2, color="grey15", max.iter = 1e6, box.padding = 0.45, max.overlaps = 100) +
  scale_color_manual(values = autumn.maple) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  coord_fixed() +
  labs(x = "portion of parsimony informative sites", y = "\npolytomy")

nCID.vs.polytomy.plot.txt <- with(tree.comp50, correlation.report(polytomy, nCID))
nCID.vs.polytomy.plot <- tree.comp50 %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  mutate(name = ifelse(gene %in% genes.to.label,gene,NA) ) %>%
  ggplot( aes(x = polytomy, y = nCID, color = process, label = name) ) +
  theme_bw() +
  theme(#legend.justification=c(0,0), legend.position=c(0.01,0.01),
        legend.position="bottom",
        axis.ticks = element_blank() ) +
  geom_smooth(method=lm, se=FALSE, color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  annotate(geom="text", x=7/8, y= 0.25/2, label = nCID.vs.polytomy.plot.txt) +
  geom_text_repel(size=2, color="grey15", max.iter = 1e6, box.padding = 0.45, max.overlaps = 100) +
  scale_color_manual(values = autumn.maple) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  coord_fixed() +
  labs(x = "polytomy", y = "tree distance\n(normalized clustering information distance)")

nCID.vs.info.plot.txt <- with(tree.comp50, correlation.report(info.sites, nCID))
nCID.vs.info.plot <- tree.comp50 %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  mutate(name = ifelse(gene %in% genes.to.label,gene,NA) ) %>%
  ggplot( aes(x = info.sites, y = nCID, color = process, label = name) ) +
  theme_bw() +
  theme(#legend.justification=c(0,0), legend.position=c(0.01,0.01),
        legend.position="bottom",
        axis.ticks = element_blank() ) +
  # geom_smooth(method=lm, se=FALSE, color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  annotate(geom="text", x=0.25/2, y= 0.25/2, label = nCID.vs.info.plot.txt) +
  geom_text_repel(size=2, color="grey15", max.iter = 1e6, box.padding = 0.45, max.overlaps = 100) +
  scale_color_manual(values = autumn.maple) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  coord_fixed() +
  labs(x = "portion of parsimony informative sites", y = "tree distance\n(normalized clustering information distance)")

misc.plots <- ggarrange(
  infoSites.vs.identity.plot, polytomy.vs.info.plot,
  nCID.vs.polytomy.plot, nCID.vs.info.plot,
  nrow = 2, ncol = 2, labels = LETTERS,
  common.legend = TRUE, legend="bottom" )

ggsave("FigS5.genetree.comp.misc.plots.pdf", misc.plots,
       width = 8, height = 8, scale = 1)

##########################
### Bootstrap cutoff of 20
##########################
tree.comp20 <- compare.gene.trees.to.organismal.tree(support.cutoff = 20)

# Combine results
tree.comp20$seq.identity <- NA
tree.comp20$info.sites <- NA
for (i in 1:(dim(tree.comp20)[1])) {
  gene.i <- tree.comp20$gene[i]
  tree.comp20$seq.identity[i] <- seq.identity[gene.i]
  tree.comp20$info.sites[i] <- info.sites[gene.i]
}

# Write the results to a file
(s <- paste0("tree.comparisons.bs20.",format(Sys.time(), "%Y.%m.%d"),".csv"))
write.csv(tree.comp20, s, row.names = FALSE)

# file.choose()
# tree.comp20 <- read.csv("tree.comparisons.bs20.2021.11.11.csv")
# tree.comp20 <- tree.comp20 %>% mutate(process=fct_relevel(tree.comp20$process, "reference"))

nCIDbyID.vs.info.plot20 <- tree.comp20 %>%
  mutate(nCID.by.ID = nCID/seq.identity ) %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  mutate(name = ifelse(gene %in% genes.to.label,gene,NA) ) %>%
  ggplot( aes(x = info.sites, y = nCID.by.ID, color = process, label = name) ) +
  theme_bw() +
  theme(legend.justification=c(0,1), legend.position=c(0.01,0.99),
        axis.ticks = element_blank() ) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(size=3, color="grey15", max.iter = 1e6, box.padding = 0.45, max.overlaps = 100) +
  scale_color_manual(values = autumn.maple) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,3.5), breaks = seq(0,3.5,0.5), labels = c("0","0.50","1.00","1.50","2.00","2.50","3.00","3.50")) +
  labs(x = "portion of parsimony informative sites", y = "\ntree distance / portion of identical sites")

polytomy.vs.info.plot20.txt <- with(tree.comp20, correlation.report(info.sites, polytomy))
polytomy.vs.info.plot20 <- tree.comp20 %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  mutate(name = ifelse(gene %in% genes.to.label,gene,NA) ) %>%
  ggplot( aes(x = info.sites, y = polytomy, color = process, label = name) ) +
  theme_bw() +
  theme(#legend.justification=c(0,0), legend.position=c(0.01,0.01),
    axis.ticks = element_blank() ) +
  geom_smooth(method=lm, se=FALSE, color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  annotate(geom="text", x=7/8, y= 7/8, label = polytomy.vs.info.plot20.txt) +
  geom_text_repel(size=2, color="grey15", max.iter = 1e6, box.padding = 0.45, max.overlaps = 100) +
  scale_color_manual(values = autumn.maple) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  coord_fixed() +
  labs(x = "portion of parsimony informative sites", y = "\npolytomy")

nCID.vs.polytomy.plot20.txt <- with(tree.comp20, correlation.report(polytomy, nCID))
nCID.vs.polytomy.plot20 <- tree.comp20 %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  mutate(name = ifelse(gene %in% genes.to.label,gene,NA) ) %>%
  ggplot( aes(x = polytomy, y = nCID, color = process, label = name) ) +
  theme_bw() +
  theme(#legend.justification=c(0,0), legend.position=c(0.01,0.01),
    legend.position="bottom",
    axis.ticks = element_blank() ) +
  geom_smooth(method=lm, se=FALSE, color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  annotate(geom="text", x=7/8, y= 1/8, label = nCID.vs.polytomy.plot20.txt) +
  geom_text_repel(size=2, color="grey15", max.iter = 1e6, box.padding = 0.45, max.overlaps = 100) +
  scale_color_manual(values = autumn.maple) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  coord_fixed() +
  labs(x = "polytomy", y = "tree distance\n(normalized clustering information distance)")

nCID.vs.info.plot20.txt <- with(tree.comp20, correlation.report(info.sites, nCID))
nCID.vs.info.plot20 <- tree.comp20 %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  mutate(name = ifelse(gene %in% genes.to.label,gene,NA) ) %>%
  ggplot( aes(x = info.sites, y = nCID, color = process, label = name) ) +
  theme_bw() +
  theme(#legend.justification=c(0,0), legend.position=c(0.01,0.01),
    legend.position="bottom",
    axis.ticks = element_blank() ) +
  geom_smooth(method=lm, se=FALSE, color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  annotate(geom="text", x=1/8, y= 1/8, label = nCID.vs.info.plot20.txt) +
  geom_text_repel(size=2, color="grey15", max.iter = 1e6, box.padding = 0.45, max.overlaps = 100) +
  scale_color_manual(values = autumn.maple) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  coord_fixed() +
  labs(x = "portion of parsimony informative sites", y = "tree distance\n(normalized clustering information distance)")

miscBS20.plots <- ggarrange(
  nCIDbyID.vs.info.plot20, polytomy.vs.info.plot20,
  nCID.vs.polytomy.plot20, nCID.vs.info.plot20,
  nrow = 2, ncol = 2, labels = LETTERS,
  common.legend = TRUE, legend="bottom" )

ggsave("FigS6.genetree.comp20.plots.pdf", miscBS20.plots,
       width = 8, height = 8, scale = 1)

##########################
### Bootstrap cutoff of 80
##########################
tree.comp80 <- compare.gene.trees.to.organismal.tree(support.cutoff = 80)

# Combine results
tree.comp80$seq.identity <- NA
tree.comp80$info.sites <- NA
for (i in 1:(dim(tree.comp80)[1])) {
  gene.i <- tree.comp80$gene[i]
  tree.comp80$seq.identity[i] <- seq.identity[gene.i]
  tree.comp80$info.sites[i] <- info.sites[gene.i]
}

# Write the results to a file
(s <- paste0("tree.comparisons.bs80.",format(Sys.time(), "%Y.%m.%d"),".csv"))
write.csv(tree.comp80, s, row.names = FALSE)

# file.choose()
# tree.comp80 <- read.csv("tree.comparisons.bs80.2021.11.11.csv")
# tree.comp80 <- tree.comp80 %>% mutate(process=fct_relevel(tree.comp80$process, "reference"))

nCIDbyID.vs.info.plot80 <- tree.comp80 %>%
  mutate(nCID.by.ID = nCID/seq.identity ) %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  mutate(name = ifelse(gene %in% genes.to.label,gene,NA) ) %>%
  ggplot( aes(x = info.sites, y = nCID.by.ID, color = process, label = name) ) +
  theme_bw() +
  theme(legend.justification=c(0,1), legend.position=c(0.01,0.99),
        axis.ticks = element_blank() ) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(size=3, color="grey15", max.iter = 1e6, box.padding = 0.45, max.overlaps = 100) +
  scale_color_manual(values = autumn.maple) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,5), breaks = seq(0,5,1), labels = c("0","1.00","2.00","3.00","4.00","5.00")) +
  labs(x = "portion of parsimony informative sites", y = "\ntree distance / portion of identical sites")

polytomy.vs.info.plot80.txt <- with(tree.comp80, correlation.report(info.sites, polytomy))
polytomy.vs.info.plot80 <- tree.comp80 %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  mutate(name = ifelse(gene %in% genes.to.label,gene,NA) ) %>%
  ggplot( aes(x = info.sites, y = polytomy, color = process, label = name) ) +
  theme_bw() +
  theme(#legend.justification=c(0,0), legend.position=c(0.01,0.01),
    axis.ticks = element_blank() ) +
  geom_smooth(method=lm, se=FALSE, color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  annotate(geom="text", x=0.25/2, y= 0.25/2, label = polytomy.vs.info.plot80.txt) +
  geom_text_repel(size=2, color="grey15", max.iter = 1e6, box.padding = 0.45, max.overlaps = 100) +
  scale_color_manual(values = autumn.maple) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  coord_fixed() +
  labs(x = "portion of parsimony informative sites", y = "\npolytomy")

nCID.vs.polytomy.plot80.txt <- with(tree.comp80, correlation.report(polytomy, nCID))
nCID.vs.polytomy.plot80 <- tree.comp80 %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  mutate(name = ifelse(gene %in% genes.to.label,gene,NA) ) %>%
  ggplot( aes(x = polytomy, y = nCID, color = process, label = name) ) +
  theme_bw() +
  theme(#legend.justification=c(0,0), legend.position=c(0.01,0.01),
    legend.position="bottom",
    axis.ticks = element_blank() ) +
  geom_smooth(method=lm, se=FALSE, color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  annotate(geom="text", x=1/8, y= 7/8, label = nCID.vs.polytomy.plot80.txt) +
  geom_text_repel(size=2, color="grey15", max.iter = 1e6, box.padding = 0.45, max.overlaps = 100) +
  scale_color_manual(values = autumn.maple) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  coord_fixed() +
  labs(x = "polytomy", y = "tree distance\n(normalized clustering information distance)")

nCID.vs.info.plot80.txt <- with(tree.comp80, correlation.report(info.sites, nCID))
nCID.vs.info.plot80 <- tree.comp80 %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  mutate(name = ifelse(gene %in% genes.to.label,gene,NA) ) %>%
  ggplot( aes(x = info.sites, y = nCID, color = process, label = name) ) +
  theme_bw() +
  theme(#legend.justification=c(0,0), legend.position=c(0.01,0.01),
    legend.position="bottom",
    axis.ticks = element_blank() ) +
  geom_smooth(method=lm, se=FALSE, color = "grey50") +
  geom_point(size = 2, alpha = 0.8) +
  annotate(geom="text", x=0.25/2, y= 0.25/2, label = nCID.vs.info.plot80.txt) +
  geom_text_repel(size=2, color="grey15", max.iter = 1e6, box.padding = 0.45, max.overlaps = 100) +
  scale_color_manual(values = autumn.maple) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  coord_fixed() +
  labs(x = "portion of parsimony informative sites", y = "tree distance\n(normalized clustering information distance)")

miscBS80.plots <- ggarrange(
  nCIDbyID.vs.info.plot80, polytomy.vs.info.plot80,
  nCID.vs.polytomy.plot80, nCID.vs.info.plot80,
  nrow = 2, ncol = 2, labels = LETTERS,
  common.legend = TRUE, legend="bottom" )

ggsave("FigS7.genetree.comp80.plots.pdf", miscBS80.plots,
       width = 8, height = 8, scale = 1)

##########################
### Stats
##########################

wilcox.test(nCID/seq.identity ~ group, data = tree.comp20)
wilcox.test(nCID/seq.identity ~ group, data = tree.comp50)
wilcox.test(nCID/seq.identity ~ group, data = tree.comp80)
# results for different bootstrap cut-offs
# BS: 20  --  W = 89, p-value = 1.004e-06
# BS: 50  --  W = 83, p-value = 5.05e-07
# BS: 80  --  W = 57, p-value = 1.645e-08

tree.comp50 %>%
  group_by(process) %>%
  summarise(n=length(gene))
# process          n
# reference       43
# other DMRT       3
# dosage comp      5
# germline sex     3
# somatic sex     17

##########################
### Taxa Info
##########################

# For each analysis taxon, how many gene trees is it represented in?
organismal.tree <- read.tree("organismal.relationships.tre")
condensed.tree.folder <- "./condensed.trees"

x <- c("id","all.trees","sex.det.genes","ref.genes")
taxon.representation <- data.frame(matrix(ncol = length(x),
                                          nrow = length(organismal.tree$tip.label)))
colnames(taxon.representation) <- x

sex.det.gene.names <- tree.comp50$gene[which(grepl("sex",tree.comp50$process) & !grepl("_",tree.comp50$gene))]
reference.gene.names <- tree.comp50$gene[grep("reference",tree.comp50$process)]

for (i in 1:length(organismal.tree$tip.label)) {
  s <- paste0("grep -l ",organismal.tree$tip.label[i]," ",condensed.tree.folder,"/*bs50.tre | grep _ -v | wc -l")
  # For example: grep -l Ephydroidea *bs50.tre | grep _ -v | wc -l
  x <- as.numeric(system(command = s, intern = TRUE))
  taxon.representation[i,] <- c(organismal.tree$tip.label[i],x,0,0)
  sum.i <- 0
  for (j in 1:length(sex.det.gene.names)) {
    s <- paste0("grep -l ",organismal.tree$tip.label[i]," ",condensed.tree.folder,"/",sex.det.gene.names[j],"*bs50.tre | grep _ -v | wc -l")
    # For example: grep -l Ephydroidea Sxl*bs50.tre | grep _ -v | wc -l
    sum.i <- sum.i + as.numeric(system(command = s, intern = TRUE))
  }
  taxon.representation[i,"sex.det.genes"] <- sum.i
  sum.i <- 0
  for (j in 1:length(reference.gene.names)) {
    s <- paste0("grep -l ",organismal.tree$tip.label[i]," ",condensed.tree.folder,"/",reference.gene.names[j],"*bs50.tre | grep _ -v | wc -l")
    # For example: grep -l Ephydroidea AATS*bs50.tre | grep _ -v | wc -l
    sum.i <- sum.i + as.numeric(system(command = s, intern = TRUE))
  }
  taxon.representation[i,"ref.genes"] <- sum.i
}

taxon.representation.plot <-
  ggtree(tr = organismal.tree) + geom_tiplab(offset = 5, hjust=0) +
  theme_tree2() +
  geom_facet(panel = "All genes",
             data = taxon.representation, geom = ggstance::geom_barh,
             aes(x = rev(as.numeric(all.trees))), fill = "gray15" , stat = "identity") +
  geom_facet(panel = "Reference genes",
             data = taxon.representation, geom = ggstance::geom_barh,
             aes(x = rev(as.numeric(ref.genes))), fill = autumn.maple[1], stat = "identity") +
  geom_facet(panel = "Sex determination genes",
             data = taxon.representation, geom = ggstance::geom_barh,
             aes(x = rev(as.numeric(sex.det.genes))), fill = autumn.maple[5], stat = "identity") +
  scale_x_continuous(limits=c(0,max(as.numeric(unlist(taxon.representation[,-1])))))

# Check the internal node numbers to find the outgroup/root number for rotation
# ggtree(tr = organismal.tree) + geom_tiplab(offset = 5, hjust=0) + geom_text(aes(label=node), hjust=-.3)
taxon.representation.plot <- ggtree::rotate(taxon.representation.plot, 58)

taxon.representation.plot

ggsave("FigS3.taxon.representation.plot.pdf", taxon.representation.plot,
       width = 6.5, height = 9, scale = 1)

##########################
### Graphical Abstract
##########################

abstract.plot <- tree.comp50 %>%
  filter(process=="somatic sex", !grepl("_",gene)) %>%
  select(gene, seq.identity, nCID ) %>%
  mutate(nCID.by.ID = nCID/seq.identity ) %>%
  arrange(desc(nCID.by.ID)) %>%
  mutate(gene = factor(gene, levels=gene)) %>%
  ggplot( aes(x = gene, y = nCID.by.ID) ) +
  theme_bw() +
  theme(legend.position="none",
        axis.ticks = element_blank(),
        panel.grid.major = element_line(colour = "#EBEBEB", size = 0.25),
        panel.grid.minor = element_line(colour = "#EBEBEB", size = 0.25),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank(),
        ) +
  geom_bar(stat="identity", fill = autumn.maple[5]) +
  labs(x = NULL, y = "tree distance / sequence identity") +
  ylim(c(0,4.5))+ coord_flip()

abstract.plot

ggsave("graphical.abstract.preliminary.pdf",abstract.plot,
       width = 2.75, height = 2.25, scale = 1.1)

##########################
### Misc. Genes
##########################
misc.comp <- compare.gene.trees.to.organismal.tree(
  input.tree.folder = "./misc.trees/",
  condensed.trees.output.folder = "./condensed.misc.trees",
  support.cutoff = 50 )

# Alignment metrics
{
  phy.files <- list.files("./misc.phy/")
  phy.files <- phy.files[grepl("\\.phy$",phy.files)]

  # Get gene names from the phylip files
  phy.names <- sub("blastp\\.\\w+\\.","",phy.files)
  phy.names <- str_split_fixed(phy.names,"\\.",2)[,1]

  seq.identity <- vector()
  info.sites <- vector()
  for (i in 1:length(phy.files)) {
    x <- read.phylip(paste0("misc.phy/",phy.files[i]))
    seq.identity[i] <- mean.pairwise.identity(x)
    info.sites[i] <- parsimony.informative.sites(x)
    cat(i,phy.names[i],"\tseq identity:",seq.identity[i],"\tinformative sites:",info.sites[i],"\n")
  }
  names(seq.identity) <- phy.names
  names(info.sites) <- phy.names
}

# Combine results
misc.comp$seq.identity <- NA
misc.comp$info.sites <- NA
for (i in 1:(dim(misc.comp)[1])) {
  gene.i <- misc.comp$gene[i]
  misc.comp$seq.identity[i] <- seq.identity[gene.i]
  misc.comp$info.sites[i] <- info.sites[gene.i]
}

# Write the results to a file
(s <- paste0("misc.comparisons.bs50.",format(Sys.time(), "%Y.%m.%d"),".csv"))
write.csv(misc.comp, s, row.names = FALSE)

# file.choose()
# misc.comp <- read.csv("misc.comparisons.bs50.2021.11.11.csv")
# misc.comp <- misc.comp %>% mutate(process=fct_relevel(misc.comp$process, "reference"))

df <- misc.comp %>%
  select(gene, process, species, taxa, polytomy, paraphyly, seq.identity, info.sites, nCID, nCID.ci ) %>%
  mutate(nCID.by.ID = nCID/seq.identity ) %>%
  arrange(desc(gene)) %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  mutate(gene = factor(gene, levels=gene)) %>%
  mutate(sample.size = paste0("n = ",taxa,"\n(",species," sp.)") ) %>%
  pivot_longer(cols = c(5:9,11))

df$nCID.ci[which(df$name != "nCID")] <- NA
df$sample.size[which(df$name != "seq.identity")] <- NA

df <- df %>%
  mutate(name = recode(as.character(name), seq.identity = "sequence\nidentity")) %>%
  mutate(name = recode(as.character(name), info.sites = "informative\nsites")) %>%
  mutate(name = recode(as.character(name), nCID = "tree distance\n(nCID)")) %>%
  mutate(name = recode(as.character(name), nCID.by.ID = "nCID/identity")) %>%
  mutate(name = factor(name, levels=unique(name)[c(3,4,1,2,5,6)]))

gene.tree.summary.stats.misc.plot <- df %>%
  ggplot( aes(x = gene, y = value, fill = process) ) +
  theme_bw() +
  theme(legend.position="bottom",
        axis.ticks = element_blank(),
        axis.text = element_text(size = 6),
        strip.text = element_text(face="bold", size = 6),
        strip.background = element_blank()) +
  facet_grid(.~name, scales = "free_x") +
  geom_bar(stat="identity") +
  scale_fill_manual(values = autumn.maple[4:5]) +
  guides(fill = guide_legend(reverse=TRUE)) +
  geom_errorbar(aes(ymin=nCID.ci, ymax=nCID.ci),
                color = "grey40", size = 1, width = 0.7) +
  geom_text(aes(x = gene, y = 1, label = sample.size),
            color = "grey65", size = 2, hjust = 1) +
  scale_y_continuous(n.breaks = 3) +
  labs(x = NULL, y = NULL) +
  coord_flip()

nCIDbyID.vs.info.misc.plot <- misc.comp %>%
  mutate(nCID.by.ID = nCID/seq.identity ) %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  ggplot( aes(x = info.sites, y = nCID.by.ID, color = process, label = gene) ) +
  theme_bw() +
  theme(legend.justification=c(0,1), legend.position=c(0.01,0.99),
        axis.ticks = element_blank() ) +
  geom_point(size = 3, alpha = 0.85) +
  geom_text_repel(size=3, color="grey15") +
  scale_color_manual(values = autumn.maple[4:5]) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,5), breaks = seq(0,5,1), labels = c("0","1.00","2.00","3.00","4.00","5.00")) +
  labs(x = "portion of parsimony informative sites", y = "\ntree distance / portion of identical sites")

info.vs.identity.misc.plot <- misc.comp %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  ggplot( aes(x = seq.identity, y = info.sites, color = process, label = gene) ) +
  theme_bw() +
  theme(#legend.justification=c(0,0), legend.position=c(0.01,0.01),
    axis.ticks = element_blank() ) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(size=2, color="grey15") +
  scale_color_manual(values = autumn.maple[4:5]) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  coord_fixed() +
  labs(x = "portion of identical sites", y = "\nportion of parsimony informative sites")

polytomy.vs.info.misc.plot <- misc.comp %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  ggplot( aes(x = info.sites, y = polytomy, color = process, label = gene) ) +
  theme_bw() +
  theme(#legend.justification=c(0,0), legend.position=c(0.01,0.01),
    axis.ticks = element_blank() ) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(size=2, color="grey15") +
  scale_color_manual(values = autumn.maple[4:5]) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  coord_fixed() +
  labs(x = "portion of parsimony informative sites", y = "\npolytomy")

nCID.vs.polytomy.misc.plot <- misc.comp %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  ggplot( aes(x = polytomy, y = nCID, color = process, label = gene) ) +
  theme_bw() +
  theme(#legend.justification=c(0,0), legend.position=c(0.01,0.01),
    legend.position="bottom",
    axis.ticks = element_blank() ) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(size=2, color="grey15") +
  scale_color_manual(values = autumn.maple[4:5]) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  coord_fixed() +
  labs(x = "polytomy", y = "tree distance\n(normalized clustering information distance)")

nCID.vs.info.misc.plot <- misc.comp %>%
  mutate(gene = ifelse(grepl("_",gene),paste0(sub("_"," (",gene),")"),gene)) %>%
  ggplot( aes(x = info.sites, y = nCID, color = process, label = gene) ) +
  theme_bw() +
  theme(#legend.justification=c(0,0), legend.position=c(0.01,0.01),
    legend.position="bottom",
    axis.ticks = element_blank() ) +
  geom_point(size = 2, alpha = 0.8) +
  geom_text_repel(size=2, color="grey15") +
  scale_color_manual(values = autumn.maple[4:5]) +
  guides(color = guide_legend(reverse=TRUE)) +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25), labels = c("0","0.25","0.50","0.75","1")) +
  coord_fixed() +
  labs(x = "portion of parsimony informative sites", y = "tree distance\n(normalized clustering information distance)")

misc.comp.plots <- ggarrange(
  gene.tree.summary.stats.misc.plot, nCIDbyID.vs.info.misc.plot,
  info.vs.identity.misc.plot, polytomy.vs.info.misc.plot,
  nCID.vs.polytomy.misc.plot, nCID.vs.info.misc.plot,
  nrow = 3, ncol = 2, labels = LETTERS,
  common.legend = TRUE, legend="bottom" )

ggsave("FigS8.genetree.misc.comp.plots.pdf", misc.comp.plots,
       width = 8, height = 12, scale = 1)


# End of script

