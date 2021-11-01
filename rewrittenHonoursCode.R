

create.folders <- function(project_paths) {


  lapply(project_paths, function (x) dir.create(x, showWarnings = FALSE))
  dir.create(
    paste0(project_paths[3],"/MDS"),
    showWarnings = F
  )

  keyfile = read_csv(project_paths[4])

  return(keyfile)
}

file_paths <- function(project.folder,analysis) {
  paths <- vector(mode = "character", length = 0)

  paths <- c(paths, paste0(
    project.folder,
    "/",
    analysis
  ))

  timeStamp <- format(
    Sys.time(),
    "%Y-%m-%d-%H%M"
  )

  paths <- c(paths, paste0(
    paths[1],
    "/featureCounts"
  ))

  paths <- c(paths, paste0(
    paths[1],
    "/DE_analysis_",
    timeStamp
  ))

  paths <- c(paths, paste0(
    paths[1],
    "/keyfile.csv"
  ))

  paths <- c(paths, paste0(paths[1], "/Mapping_stats"))


  return(paths)
}

## Summary of reads assigned to features by featureCounts (after mapping)

install_libraries <- function(libraries, installer) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  lapply(libraries, function(x) {
    if (!require(x, character.only = TRUE)) {do.call(installer, list(x))}
    library(x, character.only = TRUE)
  })
}
assignment.summary <- function(paths,keyfile) {

  #paths = project_paths
  #keyfile = keyfile

  files       <- dir(
    paste0(paths[2], "/"),
    pattern = "*.counts.summary"
    )
  count.files <- paste0(paths[2], "/", files)
  FCS         <- do.call(
    "cbind",
    lapply(
      count.files,
      function(fn) data.frame(read.delim(fn, row.names = 1))
      )
    )


  FCS            <- data.frame(t(FCS))
  row.names(FCS) <- as.data.frame(keyfile)[,1]
  FCS.t          <- t(FCS)
  FCS.plot       <- melt(FCS.t)

  prop           <- prop.table(FCS.t, margin = 2)
  prop.t         <- t(prop)
  prop.plot      <- melt(prop.t)

  plot.data      <- prop.plot

  plot_colours   <- rev(c("black", "grey", viridis_pal(option = "C")(n=12)))

  pdf(paste0(paths[5], "/percent_mapped_plot.pdf"), w = 8, h = 4)
  print(ggplot(plot.data, aes(x=Var1, y=rev(value), fill=rev(Var2))) +
    geom_bar(stat="identity") +
    xlab("Sample") +
    ylab("Percentage of reads") +
    scale_fill_manual(values = plot_colours, name = "Mapping cat") +
    text_size_theme_8_labels +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()

  ###Assigned reads as total read number

  pdf(paste0(paths[5], "/total_mapped_plot.pdf"), w = 6, h = 4)
  plot.data <- FCS.plot
  print(ggplot(plot.data, aes(x=Var2, y=value, fill=Var1)) +
    geom_bar(stat="identity") +
    xlab("Sample") +
    ylab("Number of reads") +
    scale_fill_manual(values = plot_colours, name = "Mapping cat") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
  dev.off()
}

read.data <- function(keyfile,group,paths) {
  ## Read in data
  keyfile <- as.data.frame(keyfile)
  samples <- keyfile[,1]
  sample_groups <- keyfile[, which(colnames(keyfile) == group)]

  count.files <- paste0(paths[2],"/", samples, ".counts")

  tmp <- read.delim(count.files[[1]], header=T, skip=1)
  gene.lengths <- as.numeric(tmp[,"Length"])
  names(gene.lengths) <- tmp[,"Geneid"]
  rm(tmp)

  dge <- readDGE(count.files,skip=1,
                 columns=c(1,7),
                 group=sample_groups,
                 labels=as.character(samples)
  )
  dge.raw <- dge
  gene.names <- as.character(rownames(dge$counts))
  gene.names.raw <- as.character(rownames(dge.raw$counts))

  return(dge)
}

filter.wrapper <- function (keyfile,group,paths) {
  dge        <- read.data(keyfile,group,paths)
  gene.names <- as.character(rownames(dge$counts))

  dge.raw        <- dge
  gene.names.raw <- as.character(rownames(dge.raw$counts))
  FCS_folder <- paste0(paths[3], "/GEO_featureCountsSummaries")
  dir.create(FCS_folder, showWarnings = F)

  raw.counts <- data.frame(dge$counts)
  raw.counts$AGI <- row.names(raw.counts)
  write.csv(raw.counts, paste0(FCS_folder,"/raw_counts_matrix.csv"), row.names=F)

  old.dge <- dge
  old.gene.names <- gene.names

  dge <- filtering_step(raw.counts, old.dge)
  gene.names <- as.character(rownames(dge$counts))

  return(dge)
}

filtering_step <- function(raw.counts,dge) {

  boxplot(colSums(raw.counts[,1:12]))

  #42604 genes have at least 1 count
  print("Table of genes with at least one count:")
  print(table(rowSums(raw.counts[,1:12])>0))

  #39728 genes have at least 1 count in 2 samples
  print("Table of genes with at least 1 count in 2 samples")
  print(table(rowSums((raw.counts[,1:12])>0)>2))

  #histogram of number of samples having 1 or more counts per gene
  hist(rowSums((raw.counts[,1:12])>0), breaks = 100)

  #histogram of reads per gene
  hist(log10(rowSums(raw.counts[,1:12])), breaks = 100)

  min.reads <- readline("Set number of minimum Count per Million: ")
  min.samples <- readline("Set number of minimum samples with minimum reads: ")

  min.reads <- as.numeric(min.reads)
  min.samples <- as.numeric(min.samples)

  loci.2.keep <- rowSums(cpm(dge)>min.reads)>min.samples

  print("Table of genes that passed or failed the filter")
  print(table(loci.2.keep))
  hist(
    log10(
      rowSums(
        raw.counts[
          which(rownames(raw.counts) %in%
                  names(which(loci.2.keep == TRUE))),
          1:12
          ]
        )
      ),
    breaks = 100
    )

  confirm <- readline("Accept these filters? y/n ")
  return(continue.filter(confirm, raw.counts,dge, loci.2.keep))

}

continue.filter <- function (confirm, raw.counts,dge,loci.2.keep) {
  if (confirm == "y") {
    n.tags <- sum(loci.2.keep)
    old.dge <- dge
    dge <- old.dge[loci.2.keep,]
    dge$samples$lib.size <- colSums(dge$counts)

    #make sure gene.names get updated

    return(dge)
  } else if (confirm == "n") {
    filtering_step(raw.counts,dge)
  } else {
    new.confirm <- readline("Invalid response. Accept these filters? y/n ")
    continue.filter(new.confirm,raw.counts,dge)
  }
}
