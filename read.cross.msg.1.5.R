read.cross.msg <- function(file.csv,ancfile.par2, phenofile, ancfile.par1=NULL) {
   ## Read MSG ancestry probabilities from file and return an R/qtl
   ## "cross" object.

   ## MSG ancestry probabilities are contained in a simple tabular
   ## format with individual names in first column, and labels of the
   ## form "contig:position" in the first row
   
   ##file.csv is r-qtl style csv file, which is produced by pull_thin
   ##phenofile is r-qtl style csvs phenotype file. Example provided.

##USAGE
##bc e.g.
##x <- read.cross.msg("ancestry-probs-par2.tsv.sorted.csv","ancestry-probs-par2.tsv.sorted.pulled.converted.thinned","san.pheno.csv")
##sx<-scanone(x,method="hk")
##f2 e.g.
##x <- read.cross.msg("ancestry-probs-par2.tsv.sorted.csv","ancestry-probs-par2.tsv.sorted.pulled.converted.thinned.f2_rqtl","DE14_22spheno.csv.sorted",ancfile.par1="ancestry-probs-par1.tsv.sorted.pulled.converted.thinned.f2_rqtl")
##sx<-scanone(x,method="hk")


   require(qtl)
   
   if(is.null(ancfile.par1)) backcross <- TRUE
   else backcross <- FALSE##DLS added to indicate that backcross is false if 


   ## Read MSG probabilities
   prob.AA <- as.matrix(read.table(ancfile.par2, header=TRUE, row.names=1,
                                   as.is=TRUE, check.names=FALSE))
  
   if(backcross) { 
   	   genotypes <- c("BB","BA","AA") ; alleles <- c("A","B")#Note that BB reflects backcross to parent 2 in our cross. We treat as "AA" here.
   	      	   
   	}
   else   {
       genotypes <- c("AA","AB","BB") ; alleles <- c("A","B")
       prob.BB <- as.matrix(read.table(ancfile.par1, header=TRUE, row.names=1,
                                   as.is=TRUE, check.names=FALSE))
                                   stopifnot(dim(prob.BB) == dim(prob.AA))
       } 
   #print(genotypes)
  ##Read in hard genotype calls & phenotype file

  cross <- read.cross(format="csvs", genfile=file.csv, phefile = phenofile,genotypes=genotypes,alleles =alleles,
				estimate.map = FALSE
					  )
  ##Run calc.genoprob to set up the data structure correctly. These probs will be replaced with our probs later.
  
  cross <- calc.genoprob(cross)

   ## Construct arrays of probabilities, split up by contigs
   ## (Could use abind package rather than bindfuns below)

   ## Construct data frame corresponding to R/qtl csv format
     indivs <- rownames(prob.AA)
     markers <- colnames(prob.AA)
     col.labels <- strsplit(markers, ":")
     contigs <- sapply(col.labels, "[", 1)
     contigs.f <- factor(contigs)
     bp <- as.integer(sapply(col.labels, "[", 2))
     split.bp<-split(bp,f=contigs.f)
     split.markers<-split(markers,f=contigs.f)

   split.probs <- function(p) {
       p <- split(as.data.frame(t(p)), f=contigs.f)
       lapply(lapply(p, as.matrix), t)
   }
   prob.AA <- split.probs(prob.AA)
   if(backcross) {
       prob.AB <- lapply(prob.AA, function(p) 1-p)
       bindfun <- function(pAA, pAB, markers) {
           p <- array(dim=c(dim(pAA), 2),
                      dimnames=list(indivs, markers, genotypes[1:2]))
           p[,,1] <- pAA
           p[,,2] <- pAB
           p
       }
       probs <- mapply(bindfun, prob.AA, prob.AB, split(markers, f=contigs.f), SIMPLIFY=FALSE)
   }
   else {
       prob.BB <- split.probs(prob.BB)
       prob.AB <- mapply(function(pAA, pBB) 1 - pAA - pBB, prob.AA, prob.BB, SIMPLIFY=FALSE)
       bindfun <- function(pAA, pAB, pBB, markers) {
           p <- array(dim=c(dim(pAA), 3),
                      dimnames=list(indivs, markers, genotypes))
           p[,,1] <- pAA
           p[,,2] <- pAB
           p[,,3] <- pBB
           p
       }
       #recall prob.AA is par2
       probs <- mapply(bindfun, prob.BB, prob.AB, prob.AA, split(markers, f=contigs.f), SIMPLIFY=FALSE)


   }

   ## Set the genotype probability slots, in the same way as
   ## calc.genoprob does. We may need to work out suitable values for
   ## some of those attributes.

   step <- 0 ; off.end <- 0 ; stepwidth <- "fixed"

   for(contig in levels(contigs.f)) {
       
       if(backcross) 
	   cross$geno[[contig]]$prob <- probs[[contig]]
	   	   
       else {
 	   	   #check contig, if not X, just proceed
 	       if (contig != 'X')
                   cross$geno[[contig]]$prob <- probs[[contig]]
               else
                   #grab only probs[,,1] & probs[,,3]
                   cross$geno[[contig]]$prob <- probs[[contig]][,,c(1,3)]
 	   	   
             }

       contig.map<-as.numeric(split.bp[[contig]])
       attr(contig.map,"names") <-split.markers[[contig]]
       attr(cross$geno[[contig]]$prob, "map") <- contig.map
       attr(cross$geno[[contig]]$prob, "error.prob") <- 1e-04 #fixed this prob slot in v. 1.5, now can run scantwo
       attr(cross$geno[[contig]]$prob, "step") <- step
       attr(cross$geno[[contig]]$prob, "off.end") <- off.end
       attr(cross$geno[[contig]]$prob, "map.function") <- "haldane" 
       attr(cross$geno[[contig]]$prob, "stepwidth") <- stepwidth
       
   }
   cross
}
