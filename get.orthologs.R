# A function to retrieve orthologs of a given Accession
# from taxa listed in a file with NCBI taxids.
# Based on BLASTp to a local NCBI NR BLAST database.
#
# To modify the blastp CSV output format (outfmt.string) check the NCBI docs
# https://www.ncbi.nlm.nih.gov/books/NBK279684/

get.orthologs <- function(
  query.accession,
  query.user.description,
  taxids.file ="taxids.csv",
  blastdb.path = "/storage/blastdb/2021-09-29/nr",
  evalue.cutoff = 1e-15,
  max.target.seqs = 1,
  max.hsps = 1,
  num_threads = 12,
  outfmt.string = "qacc staxids sacc evalue bitscore length nident mismatch gaps stitle",
  return.results = FALSE,
  align.sequences = TRUE
)
{ # Begin the function

  # Note the time
  start.time <- Sys.time()

  # Check packages
  required.packages <- as.character(c("stringr","rentrez","phylotools"))
  for (i in 1:length(required.packages)) {
    if (!(required.packages[i] %in% rownames(installed.packages()))) {
      stop(paste0("Package missing. First, try running `install.packages('",required.packages[i],"')`"))
    }
    require(package = required.packages[i], character.only = TRUE)
  }

  # Validate input
  if (missing(query.user.description)) {
    query.user.description <- ""
    while (nchar(query.user.description) < 1) {
      x <- rentrez::entrez_summary(db="protein", id=query.accession)
      s <- paste0("Please provide a short description and press [Enter]:\n",
                  "Avoid by setting `query.user.description`\n",
                  "Entrez title: ", x$title,"\n")
      query.user.description <- readline(prompt = s)
    }
  }
  if (!grepl("sacc",outfmt.string)) {
    stop("The arguement `outfmt.string` must include `sacc`.")
  }
  sacc.col <- grep("sacc",unlist(str_split(outfmt.string," ")))

  # Initialize
  temp.query.filename <- paste0("tmp.query.file.",query.user.description,".",format(Sys.time(), "%y%m%d.%s"),".tmp")

  # Load the file of taxon IDs
  taxid.table <- read.csv(taxids.file)
  taxids <- taxid.table[,1]
  species.names <- taxids
  if (any(grepl("species",names(taxid.table)))) {
    species.names <- taxid.table$species
    if (any(grepl("order",names(taxid.table))) & any(grepl("family",names(taxid.table)))) {
      species.names <- paste0(species.names," (",taxid.table$order,": ",taxid.table$family,")")
    }
  }

  # Write the query to a file
  write(query.accession, temp.query.filename)

  # Set-up long-format output file
  output.filename.long.format <- paste0(query.user.description,".",query.accession,".",format(Sys.time(), "%y%m%d"),".blastp.results.csv")
  s <- sub("^\\d+ ","",outfmt.string)
  s <- paste0(gsub(" ",",",s),",,,,")
  write(s, output.filename.long.format)

  # Main loop
  for (i in 1:length(taxids)) {
    cat(paste(i,"\tSearching ",species.names[i],"\tfor ",query.user.description,"...\t"))

    results <- read.csv(output.filename.long.format, header=FALSE, stringsAsFactors = FALSE, row.names = NULL)
    initial.length <- dim(results)[1]

    # System call
    s <- paste0(" -query ",temp.query.filename,
                " -db ",blastdb.path,
                " -evalue ",evalue.cutoff,
                " -max_target_seqs ",max.target.seqs,
                " -max_hsps ",max.hsps,
                " -num_threads ",num_threads,
                " -taxids ",taxids[i],
                " -outfmt '10 ",outfmt.string,"'",
                " >> ",output.filename.long.format)

    system2(command = "blastp", args = s,
            wait = TRUE, stderr = NULL)

    results <- read.csv(output.filename.long.format, header=TRUE, stringsAsFactors = FALSE, row.names = NULL)

    if (dim(results)[1] == initial.length) {
      cat("found:",results[sacc.col,initial.length],"\n")
    } else {
      cat("no hit\n")
    }

  } # End main loop

  # Delete temp file
  system2(command = "rm", args = temp.query.filename)

  # Set-up short-format output file: just the list of subject accession numbers
  output.filename.short.format <- sub("results\\.csv","accessions.txt",output.filename.long.format)
  x <- paste0(results$sacc, collapse = "\n")
  write(x, output.filename.short.format)

  # Remove any blank lines in the file
  system2(command = "sed", args = paste0("'/^$/d' -i ",output.filename.short.format))

  # Display run info
  cat(paste0("Completed BLASTp: ",query.accession," ",query.user.description,"\n"))
  cat(paste0("Searched: ",i," taxa\n"))
  x <- dim(results)[1]
  cat(paste0("Found: ",x," matches\n"))

  if (align.sequences) {
    create.alignment.from.acc.list(
      accession.list.file = output.filename.short.format,
      num_threads = num_threads
    )
    cat(paste0("Completed alignment: ",query.accession," ",query.user.description,"\n"))
  }

  cat(paste0("Run time using ",num_threads," threads: ",format(Sys.time() - start.time),"\n"))

  if (return.results) {
    return(results)
  }

} # End function


# A function to align sequences starting from a list of NCBI Accessions
# Uses Clustal-Omega and edits sequence names to conform to RAxML expectations

create.alignment.from.acc.list <- function(
  accession.list.file,
  output.filename.base = NULL,
  save.only.phy.output = FALSE,
  num_threads = 8,
  verbose = TRUE
)
{
  # Check packages
  required.packages <- as.character(c("stringr","rentrez","phylotools"))
  for (i in 1:length(required.packages)) {
    if (!(required.packages[i] %in% rownames(installed.packages()))) {
      stop(paste0("Package missing. First, try running `install.packages('",required.packages[i],"')`"))
    }
    require(package = required.packages[i], character.only = TRUE)
  }

  # Read the accession list file
  subject.ids <- unlist(read.delim(accession.list.file, sep = "\n", header = FALSE))
  if (subject.ids[1]=="sacc") { subject.ids <- subject.ids[-1] }

  if (is.null(output.filename.base)) {
    output.filename.phy <- sub("txt$","phy",accession.list.file)
    output.filename.phy <- sub("accessions","ortholog.seqs",output.filename.phy)
  } else {
    output.filename.phy <- paste0(output.filename.base,".phy")
  }
  output.filename.faa <- sub("\\.phy",".faa",output.filename.phy)
  output.filename.dnd <- sub("\\.phy",".dnd",output.filename.phy)
  output.filename.aligned.faa <- sub("\\.faa",".aligned.faa",output.filename.faa)

  if (file.exists(output.filename.phy)) {
    x <- readline(prompt = paste0(output.filename.phy,"  already exists. Overwrite? (y/n) "))
    if (x =="y") {
      system2(command = "rm", args = output.filename.phy)
    } else {
      x <- readline(prompt = "Append? (y/n) -- 'No' will quit. ")
      if (x !="y") {
        stop("Stopped.")
      }
    }
  }
  if (file.exists(output.filename.faa)) {
    x <- readline(prompt = paste0(output.filename.faa,"  already exists. Overwrite? (y/n) "))
    if (x =="y") {
      system2(command = "rm", args = output.filename.faa)
    } else {
      x <- readline(prompt = "Append? (y/n) -- 'No' will quit. ")
      if (x !="y") {
        stop("Stopped.")
      }
    }
  }

  # Work around for problematic accessions. (This is a stupid hack!)
  # Some accessions cause entrez_summary to throw an error that stops the parent function.
  # As I discover them, they're simply singled out here and provided manual species annotations.
  # Most seem to be from Portunus trituberculatus, so I've removed this species from the taxid list.
  problem.IDs <- c("MPC72813","MPC45793","RXG54398","MPC11630","MPC25974")
  names(problem.IDs) <- c("Portunus_trituberculatus","Portunus_trituberculatus","Armadillidium_vulgare","Portunus_trituberculatus")

  # Main loop
  for (i in 1:length(subject.ids)) {
    # Get sequence
    s <- rentrez::entrez_fetch(db="protein", id=subject.ids[i], rettype="fasta")
    s <- sub("\\n\\n","",s)
    seq.i <- str_split_fixed(s,"\\n",2)[1,2]

    # Get subject species name
    if (!(subject.ids[i] %in% problem.IDs)) {
      x <- rentrez::entrez_summary(db="protein", id=subject.ids[i])
      x <- rentrez::entrez_summary(db="taxonomy", id=x$taxid)
      species.i <- paste0(x$genus,"_",x$species)
    } else {
      species.i <- names(problem.IDs)[which(problem.IDs == subject.ids[i])]
    }

    # Write output
    s <- paste0(">",species.i,"_",subject.ids[i],"\n",seq.i)
    write(s, output.filename.faa, append = TRUE)

    # Report progress
    cat(i,"\tPulled sequence: ",gsub("_"," ",species.i),subject.ids[i],"\n")
  } # End main loop

  # Clustal-Omega alignment
  s <- paste0(" -i ",output.filename.faa,
              " --iter=2 --force",
              " -o ",output.filename.aligned.faa,
              " --outfmt=fa --wrap=60 --output-order=tree-order",
              " --guidetree-out=",output.filename.dnd,
              " --threads=",num_threads )
  if (verbose) {
    s <- paste0(s," --verbose")
  }
  cat("Performing Clustal-Omega alignment.\n")
  system2(command = "clustalo", args = s,
          wait = TRUE, stderr = NULL)

  # Convert aligned FASTA file to PHYLIP format
  convert.fa.to.phy(file = output.filename.aligned.faa, outfile = output.filename.phy)

  # Optional clean-up
  if (save.only.phy.output) {
    system2(command = "rm", args = output.filename.faa)
    system2(command = "rm", args = output.filename.aligned.faa)
    system2(command = "rm", args = output.filename.aligned.dnd)
  }

  cat("\n")
  # cat("\nDone.\n")

} # End Function


convert.fa.to.phy <- function(file, outfile) {
  require(phylotools)
  if (missing(outfile)) {
    outfile <- sub("\\.aligned","",file)
    outfile <- sub("\\.faa",".phy",outfile)
  }
  x <- read.fasta(file = file)
  dat2phylip(x, outfile = outfile)
}

