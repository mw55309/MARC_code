#!/usr/bin/env Rscript

library(poRe)
source("get.kmer.events.R")

args <- commandArgs(trailingOnly = TRUE)
#f5 <- "/nanopore/MARC/DD/DD_575_R7.3_MARC_K12_Ia_04_24_15_Run_2/pass/DEAMERNANOPORE_DD_575_R7.3_MARC_1a_Run2_0306_1_ch308_file9_strand.fast5"
f5 <- args[1]

# 2D Fastq
fq <- get_fastq(f5, which="2D")

# link events to 2D kmers
aln <- get.kmer.events(f5)

# get 2D sequence
four <- strsplit(fq$'2D', "\n")
seq <- four[[1]][2]

# add 2D kmer position to the alignment table
kmer.length <- 5
seq.pos     <- 1
scan        <- 10

kmer.pos    <- rep("NA", length=nrow(aln))

# iterate over aln object
for (i in 1:nrow(aln)) {

	current.kmer <- aln$kmer[i]
	
	current.pos  <- seq.pos
	kmer.2d      <- substr(seq, seq.pos, seq.pos + kmer.length - 1)
	while(kmer.2d != current.kmer && (current.pos - seq.pos) <= scan) {
		current.pos <- current.pos + 1
		kmer.2d     <- substr(seq, current.pos, current.pos + kmer.length - 1)
	}

	# if we have 
	if (kmer.2d == current.kmer) {
		seq.pos <- current.pos
		kmer.pos[i] <- seq.pos
	}
	
}

aln$kmer.pos <- kmer.pos

write.table(aln, "f5.kmer.withpos.events", row.names = FALSE, quote = FALSE, sep="\t")





