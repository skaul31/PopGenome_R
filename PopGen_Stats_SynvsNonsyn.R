library(PopGenome)

# Use wormbase to identify the physical start and stop position of the gene of interest

# Read in your gene of interest along with the gff file - using par-2 as example 

GENOME.classpar2 <- readVCF("Elegans.hf.isotype.vcf.gz", numcols=10000, tid="III", from=1081402, to=1095082, approx=FALSE, include.unknown=TRUE, out="", parallel=FALSE, gffpath="WS283.gff3")

# In order to get statistics for synonymous and nonsynonymous SNPs we need to set the snps using a reference file 

# I used the coding sequence of par-2 with the largest isoform

# make sure to set save.codons = TRUE as this allows for the synonymous and nonsynonymous SNPs to be saved as a dataframe 

GENOME.classpar2 <- set.synnonsyn(GENOME.classpar2, ref.chr="ChrIII.fasta", save.codons = TRUE)


# Now calculate neutrality stats for each set of sites 

GENOME.class_nonsyn <- neutrality.stats(GENOME.classpar2, subsites="nonsyn", FAST=TRUE)

GENOME.class_syn <- neutrality.stats(GENOME.classpar2, subsites="syn", FAST=TRUE)

# Number of segregating synonymous and nonsynonymous sites are stored in n.segregating.sites 

GENOME.class_syn@n.segregating.sites

GENOME.class_nonsyn@n.segregating.sites


# Store Tajima's D values as separate vectors for each set of sites 

synTaj <- GENOME.class_syn@Tajima.D

nonsynTaj <- GENOME.class_nonsyn@Tajima.D

synTaj

nonsynTaj

# For synonymous and nonsynymous pi similar code can be used but instead of neutrality.stats use diversity.stats 

GENOME.class_synD <- diversity.stats(GENOME.classpar2, subsites="syn")

GENOME.class_nonsynD <- diversity.stats(GENOME.classpar2, subsites="nonsyn")

# To get pi for each set of sites just access the value with nuc.diversity.within

GENOME.class_synD@nuc.diversity.within

GENOME.class_nonsynD@nuc.diversity.within



