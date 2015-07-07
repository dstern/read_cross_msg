plot.genotypes.qtl <- function(msg_cross_object, scanone_object, phenotype, format, zoom, plot_title){

####plot.genotypes.qtl plots all the genotypes of individuals ordered according to the value of a phenotype and the QTL corresponding to that phenotype	
####USAGE
####plot.genotypes.qtl(msg_cross_object, scanone_object, phenotype, format, zoom, plot_title)
#####msg_cross_object - an object generated with read.cross.msg.1.5.R code
#####scanone_object - an object generated applying scanone function to msg_cross_object
#####phenotype - the name of the phenotype column in msg_cross_object whose QTL will be plotted (character)
#####format - It can take 2 values, "genome" to plot the entire QTL map or "region" to plot a zoom of the genome (character)
#####zoom - A vector with 3 values of the form c("2", 1e7, 2e7). The first value is the chromosome name as character, the second and third values are the starting and ending zoom positions.
#####plot_title - The title for the plot

	
    ##########
    ##########
    ###To change the relative sizes of the QTL and genotype graphs among each other. This can be modified by the user
    geno_plot_size <- 0.8   #Default values
    qtl_plot_size <- 0.55
    

	###########
	###########
	###To Test COLORS
	#cols = rev(rainbow(101, s = 1, v = 1, start = 0, end = 0.3)) #Gradient of 100 rainbow colors
	cols = rev(heat.colors(101)) #Gradient of 100 heat colors, default.
	cols[1] = "white" 
	#cols[51] = "green"
	X_cols = rainbow(101, s = 1, v = 1, start = 0.675, end = 0)
	
	#par(fig = c(0, 1, 0, 0.5)) #To test colors for autosomal chromosomes
	#plot(c(0, 101), c(0,1))
	#rect(seq(1,100,1), 0, seq(2,101,1), 1, col= cols[])
	#title("Autosomal colors")
	
	#par(fig = c(0, 1, 0.5, 1), new = TRUE) #To test colors for X chromosome
	#plot(c(0, 101), c(0,1))
	#rect(seq(1,100,1), 0, seq(2,101,1), 1, col= X_cols[])
	#title("X colors")
	
	############
	############
	###GENOTYPES
	markers <- NULL 
	markers_start <- NULL
	markers_end <- NULL
	probability	<- NULL
	toPlot_probability <- NULL
	priors <- NULL

	############
	############
	###Number of CHROMOSOMES and their NAMES	
	chr_names <- rownames(summary(msg_cross_object$geno))  #Number and names of chromosomes
	#print(chr_names)
	#length(chr_names)
	
	chr_pos <- matrix(0,2,length(chr_names))  #chr_pos will contain the end and start positions of the chromosomes
	chr_pos
	chr_lengths <- NULL
	
	for(i in 1:length(chr_names)){
		
		chr_lengths[i] <- strsplit(as.character(names(msg_cross_object$geno[chr_names[i]][[1]]["map"][[1]][length(msg_cross_object$geno[chr_names[i]][[1]]["map"][[1]])])),":")[[1]][2]  #Takes the first and last marker on each chromosome
		chr_lengths <- as.numeric(chr_lengths)

		chr_pos[1,i] <- sum(chr_lengths[1:i-1]) #Makes a 2xnumber_of_chromosomes matrix with the start position in the 1st row and the end position en the 2nd row for each chromosome (columns)
		chr_pos[2,i] <- sum(chr_lengths[1:i])
	
	}
	
	#############
	#############
	###PHENOTYPES
	
	###Number and names of individuals in the experiments
	fly_names <- as.character(msg_cross_object$pheno["id"][,1])
	
	###To sort the individuals by phenotype value in increasing order (eg. by IPI)
	ord_phen <- msg_cross_object$pheno[order(msg_cross_object$pheno[,phenotype]),] #Sorts the individuals by phenotype
	temp <- which(ord_phen[,phenotype] != 0) #Takes the position of phenotypes, this is, individuals without a NAs
	ord_phen <- ord_phen[temp,] #Winnows NAs
	phen_ordered_pos <- as.numeric(rownames(ord_phen)) #Takes the positions in the original read.cross.msg object of individuals sorted by phenotype
	
	fly_names <- fly_names[phen_ordered_pos] #Array with names of individuals with phenotype sorted by phenotype
	
	##################
	##################
	#####PLOT
	
	width_chr_plot <- seq(1/length(fly_names),1,1/length(fly_names)) #To set the width of each chromosme image in the plot according to the amount of individuals
	width_chr_one <- width_chr_plot[1]


	##################
	########plot WHOLE genome
	
	if(format == "genome"){
		
		par(fig = c(0, 1, 0, geno_plot_size)) #Divides the graphical display in two, this is for the genotypes plot
		
		plot(c(0, chr_pos[2,dim(chr_pos)[2]]), c(0, 1), type = "n", xlab="Chr", ylab=phenotype, xaxt = "n", yaxt = "n")
		
		axis(1, at = c(chr_pos[1,]), labels = as.expression(chr_names)) #Draws X axis, chromosome positions
		chr_names_lengths <- chr_lengths
		names(chr_names_lengths) <- chr_names
		print("Chromosome Lengths")
		print(chr_names_lengths) 
		
		#Gets phenotype references for the Y axis
		y_ref <- c(msg_cross_object$pheno[phen_ordered_pos[1],phenotype], msg_cross_object$pheno[phen_ordered_pos[round(length(phen_ordered_pos)*0.2)],phenotype], msg_cross_object$pheno[phen_ordered_pos[round(length(phen_ordered_pos)*0.4)],phenotype], 
		msg_cross_object$pheno[phen_ordered_pos[round(length(phen_ordered_pos)*0.6)],phenotype], msg_cross_object$pheno[phen_ordered_pos[round(length(phen_ordered_pos)*0.8)],phenotype], msg_cross_object$pheno[phen_ordered_pos[round(length(phen_ordered_pos)*1)],phenotype])
		#print(y_ref)
		
		axis(2, at = c(0,0.2, 0.4, 0.6, 0.8, 1), labels = y_ref)
		
		for(j in 1:length(fly_names)){ #Loop to take each of the ordered by phenotype individuals
			for(i in 1:length(chr_names)){ #Loop to enter in all the chromosome of an individual
				
				markers	<- NULL
				markers_start <- NULL
				markers_end	<- NULL
				probability	<- NULL
				toPlot_probability <- NULL
				priors	<- NULL
				
				####Markers positions
				markers <- unlist(strsplit(names(msg_cross_object$geno[chr_names[i]][[1]]["prob"][[1]][,,"BB"][fly_names[j],]),":"))
				#print(markers)
				
				markers_end <- as.numeric(markers[seq(2, length(markers), 2)]) + chr_pos[1,i]
				#print(markers_end)
				
				markers_start[1] <- chr_pos[1,i]
				markers_start[2:length(markers_end)] <- markers_end[1:length(markers_end) - 1] 
				#markers_start
				
				####Probabilities of each marker
				probability <- msg_cross_object$geno[chr_names[i]][[1]]["prob"][[1]][,,"BB"][fly_names[j],]
				toPlot_probability <- round(probability, digits = 2) * 100 + 1
				
				####Priors
				priors <- which(probability == 0.5)
				#print("priors")
				#print(priors)
				####Plot Chromosomes
				if(chr_names[i] == "X"){    #If the marker is located on chromosome X and its probability to come from parent 1 is 0, plot it in blue
						rect(markers_start, width_chr_plot[j], markers_end, width_chr_plot[j] + width_chr_one, border = X_cols[toPlot_probability], col = X_cols[toPlot_probability])
						if(length(priors) != 0 ){   #Plot priors in black if there is any
							rect(markers_start[priors], width_chr_plot[j], markers_end[priors], width_chr_plot[j] + width_chr_one, border = "black", col = "black")
							}
					}
					
				else{   #Plot priors in black if there is any
					rect(markers_start, width_chr_plot[j], markers_end, width_chr_plot[j] + width_chr_one, border = cols[toPlot_probability], col = cols[toPlot_probability])
					if(length(priors) != 0 ){
						rect(markers_start[priors], width_chr_plot[j], markers_end[priors], width_chr_plot[j] + width_chr_one, border = "black", col = "black")
					}
				}
				
			}
		}
		################
		################
		###QTL PLOTTING
		par(fig = c(0, 1, qtl_plot_size, 1), new = TRUE) #Divides the graphical display in two, this is for the QTL map
	
		genome_pos_scanone <- NULL #Contains the position of each marker in the genome
        
        if(dim(scanone_object)[2] == 3){    #To avoid errors when there is just one phenotype in the scanone object
            phenotype_qtl = "lod" 
        }
        else{
            phenotype_qtl = phenotype
        }

		for(i in 1:length(rownames(scanone_object))){
	
			genome_pos_scanone[i] <- as.numeric(strsplit(rownames(scanone_object), ":")[[i]][2]) + chr_pos[1,which(chr_names == strsplit(rownames(scanone_object), ":")[[i]][1])] #Converts the marker positions in the chromosome to positions in the genome
						
			}

		plot(genome_pos_scanone,scanone_object[,phenotype_qtl], type = "l", xlab = "", ylab="LOD score", main = plot_title, xaxt = "n")
		axis(1, at = c(chr_pos[1,]), labels = rep("",length(chr_pos[1,]))) #Draws X axis, chromosome positions
		
	}
	
	##################
	########plot ZOOM region
	
	if(format == "region"){ ####################PLOT ZOOM
			
		chr_usr <- zoom[1]  #Chromosome
		#print(chr_usr)
		#chr_usr <- which(chr_names == chr_usr)
		#print(zoom)
		#print(chr_usr)
			
		gpos_1 <- as.numeric(zoom[2])   #Zoom regions
		gpos_2 <- as.numeric(zoom[3])
		
			
		####Markers positions for that chromosome
		markers <- unlist(strsplit(names(msg_cross_object$geno[chr_usr][[1]]["prob"][[1]][,,"BB"][fly_names[1],]),":")) 
		#print(markers)

		markers_end <- as.numeric(markers[seq(2, length(markers), 2)])
		#print(markers_end)

		markers_start[1] <- 0
		markers_start[2:length(markers_end)] <- markers_end[1:length(markers_end) - 1] 
		#print(markers_start)
			
		pos_1 <- which(abs(markers_start - gpos_1) == min(abs(markers_start - gpos_1)))  ###Lower interval
		#print(markers_start[pos_1])
		pos_2 <- which(abs(markers_end - gpos_2) == min(abs(markers_end - gpos_2)))  ###Higher interval
		#print(markers_end[pos_2])
			
		#####Take User Defined Positions
		start_user <- markers_start[pos_1]
		end_user <- markers_end[pos_2]
		#print(start_user)
		#print(end_user)
		
		par(fig = c(0, 1, 0, geno_plot_size)) #Divides the graphical display in two, this is for the genotypes plot
			
		plot(c(markers_start[pos_1], markers_end[pos_2]), c(0, 1), type = "n", xlab= paste("chr", chr_usr, sep = " "), ylab=phenotype, yaxt = "n")#, xaxt = "n", yaxt = "n") #Draws the display
		
		#Gets phenotype references for the Y axis
		y_ref <- c(msg_cross_object$pheno[phen_ordered_pos[1],phenotype], msg_cross_object$pheno[phen_ordered_pos[round(length(phen_ordered_pos)*0.2)],phenotype], msg_cross_object$pheno[phen_ordered_pos[round(length(phen_ordered_pos)*0.4)],phenotype], 
		msg_cross_object$pheno[phen_ordered_pos[round(length(phen_ordered_pos)*0.6)],phenotype], msg_cross_object$pheno[phen_ordered_pos[round(length(phen_ordered_pos)*0.8)],phenotype], msg_cross_object$pheno[phen_ordered_pos[round(length(phen_ordered_pos)*1)],phenotype])
		#print(y_ref)
		
		axis(2, at = c(0,0.2, 0.4, 0.6, 0.8, 1), labels = y_ref)
		
		########
		##ADD PLOT FEATURES

		for(j in 1:length(fly_names)){ #Loop to take each of the ordered by phenotype individuals
				
			probability	<- NULL
			toPlot_probability <- NULL
			priors	<- NULL		
				
			####Probabilities of each marker
			probability <- msg_cross_object$geno[chr_usr][[1]]["prob"][[1]][,,"BB"][fly_names[j],]
			priors <- which(probability == 0.5)
			toPlot_probability <- round(probability[pos_1:pos_2], digits = 2) * 100 + 1
				
			if(chr_usr == "X"){ #If the marker is located on chromosome X and its probability to come from parent 1 is 0, plot it in X colors
				rect(markers_start[pos_1:pos_2], width_chr_plot[j], markers_end[pos_1:pos_2], width_chr_plot[j] + width_chr_one, border = X_cols[toPlot_probability], col = X_cols[toPlot_probability])
					
				if(length(priors) > 0){ #Plot priors in black if there is any
					rect(markers_start[priors], width_chr_plot[j], markers_end[priors], width_chr_plot[j] + width_chr_one, border = "black", col = "black")
					}
				}
			
			else{
				rect(markers_start[pos_1:pos_2], width_chr_plot[j], markers_end[pos_1:pos_2], width_chr_plot[j] + width_chr_one, border = cols[toPlot_probability], col = cols[toPlot_probability])
						
				if(length(priors) > 0){ #Plot priors in black if there is any
					rect(markers_start[priors], width_chr_plot[j], markers_end[priors], width_chr_plot[j] + width_chr_one, border = "black", col = "black")
					}		
				}
				

			}
			
		#######
		#######
		#######
		##QTL PLOTTING
		
		par(fig = c(0, 1, qtl_plot_size, 1), new = TRUE) #Divides the graphical display in two, this is for the QTL map
	
		chr_pos_scanone <- which(scanone_object$chr == chr_usr)[pos_1:pos_2]
		
		scanone_object[chr_pos_scanone,]
		

        
        if(dim(scanone_object)[2] == 3){    #To avoid errors when there is just one phenotype in the scanone object
            phenotype_qtl = "lod" 
        }
        else{
        phenotype_qtl = phenotype
        }


		plot(scanone_object[chr_pos_scanone,"pos"],scanone_object[chr_pos_scanone,phenotype_qtl], type = "l", xlab = "", ylab="LOD score",xaxt = "n")
		
		plot_title2 = paste("zoom ", gpos_1, "bp - ", gpos_2, "bp ",sep = "")
		plot_title3 = paste("chr ", chr_usr, sep = "")

		title(plot_title, line = 3)
		title(plot_title2, line = 2, cex.main = 0.8)
		title(plot_title3, line = 1, cex.main = 0.8)
			
		}
		
		if(format != "genome" & format != "region"){
			print("Choose between plotting the whole genome or a user defined region")
			}
			
	
	}#END
