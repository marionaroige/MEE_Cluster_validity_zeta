####################################################
# Mariona Roigé Valiente, 2016                     #
# Ecological Informatics, BPRC, Lincoln University #              			  
# mariona.roige.valiente@gmail.com       	   #
####################################################

setwd('/home/mariona/Documents/00_PhD/0_Thesis_outline/03_Zetadiversity/02_paperMME/Materials for online')
library(zetadiv)
library(plyr)
library(dplyr)
library(stringr)

# Import data
# data are organized in a matrix (regions, species). Regions are in the row names. 
data <- read.table('dataclean2006.txt',header = TRUE, stringsAsFactors = FALSE) 
cells_raw <- read.table('cells.csv', sep=',',na='',stringsAsFactors = FALSE)
# cells are the output from sM.codebook. First column is the cell(neuron) name. 

nreg <- dim(data)[1] # number of regions in the original dataset
nspe <- dim(data)[2] # number of species in the original dataset
ncell <- dim(cells_raw)[1] # number of cells in the SOM output map
max_reg <- (dim(cells_raw)[2]-1) #(Number of columns-1) tells the maximum number of regions allocated in one cell 
empty_cells <- which(is.na(cells_raw[,2])) # Identify out of the total number of cells which ones are empty
single_cells <- which(is.na(cells_raw[,3])&!is.na(cells_raw[,2])) # Identify out of the total number of cells which have only one region clustered in them
non_empty_cells <- ncell-length(empty_cells)
cells <- cells_raw[,-1]
row.names(cells) <- cells_raw[,1]

# Create subsets of the data per cluster to compute its zeta values wih function Zeta.decline()

# 1- Check which sites are in row i in 'cells'
# 2- Extract the row from 'data' for each one of these i sites
# 3- Delete in each matrix the extra NA rows created


listofmates <- list()
matsites <- data.frame(matrix(, nrow = max_reg, ncol = nspe)) 										

for (i in 1:ncell){	 

	sites <- cells[i,]
	matsites <- data.frame(matrix(, nrow = max_reg, ncol = nspe))
	
		for (j in 1:length(sites)){
		
			sitetosearch <- sites[j]
			test <- is.na(sitetosearch)
				
			if (test == FALSE){ 
			pos <- grep(sitetosearch, row.names(data))
			tostore <- data[pos,]
			matsites[j,] <- as.matrix(tostore)
			row.names(matsites)[j] <- row.names(tostore)
					  }}
matsites <- na.omit(matsites) 
names(matsites) <- names(data)
listofmates[[i]] <- matsites}
summary(listofmates) # Summary should output a list of ncell data.frames each one of length nspe

# Eliminate from the list those cells that have none or only
# one region clustered in them and thus for which zetadiv cannot be computed. 

pos_to_delete <- sort(c(single_cells,empty_cells))
clusters_ready <- listofmates[-pos_to_delete]

##########################################################################################################################
# A version of the function Zeta.decline() to be able to run it trough lapply to datasets of diffrent size that         ##
# are sometimes smaller than 10 regions. 										##
# Changes from the original: Zeta orders are specified as the number of rows of the data, id est, the number of regions.##
# Also, plotting has been shut off.										        ##
##########################################################################################################################



Zeta.decline.v2 <- function (data.spec, orders = 1:dim(data.spec)[1], sam = 1000, plot = FALSE, 
    sd = TRUE) 
{
    if (class(data.spec) != "data.frame") {
        stop(paste(deparse(substitute(data.spec)), " is a ", 
            class(data.spec), ". It must be a data frame.", sep = ""))
    }
    if (max(orders) > dim(data.spec)[1]) {
        stop("Wrong value for \"orders\": the maximum value must be equal or lower than the number of sites.")
    }
    x <- dim(data.spec)[1]
    zeta.val <- numeric()
    zeta.val.sd <- numeric()
    for (j in orders) {
        if (choose(x, j) > sam) {
            u <- rep(NA, sam)
            for (z in 1:sam) {
                samp <- sample(1:x, j, replace = FALSE)
                u[z] <- sum(apply(data.spec[samp, ], 2, prod))
            }
        }
        else {
            u <- rep(NA, choose(x, j))
            samp <- combn(1:x, j)
            for (z in 1:dim(samp)[2]) {
                u[z] <- sum(apply(data.spec[samp[, z], ], 2, 
                  prod))
            }
        }
        zeta.val[j] <- mean(u)
        zeta.val.sd[j] <- sd(u)
    }
    zeta <- list()
    zeta$zeta.order <- orders
    zeta$zeta.val <- zeta.val
    zeta$zeta.val.sd <- zeta.val.sd
    zeta.val.log <- log10(zeta.val)
    zeta.val.log[which(is.infinite(zeta.val.log))] <- NA
    zeta.exp <- lm(zeta.val.log ~ c(orders), na.action = na.omit)
    zeta$zeta.exp <- zeta.exp
    zeta.pl <- lm(zeta.val.log ~ log10(c(orders)), na.action = na.omit)
    zeta$zeta.pl <- zeta.pl
    zeta$aic <- AIC(zeta$zeta.exp, zeta$zeta.pl)
    if (plot == TRUE) {
        par(mfrow = c(1, 3))
        if (sd == TRUE) {
            plot(orders, zeta.val, xlab = "Zeta order", ylab = "Zeta-diversity", 
                pch = 20, ylim = c(0, zeta.val[1] + zeta.val.sd[1]), 
                main = "Zeta diversity decline")
            lines(orders, zeta.val)
            for (i in orders) {
                suppressWarnings(arrows(i, zeta.val[i], i, zeta.val[i] + 
                  zeta.val.sd[i], angle = 90, length = 0.1))
                suppressWarnings(arrows(i, zeta.val[i], i, zeta.val[i] - 
                  zeta.val.sd[i], angle = 90, length = 0.1))
            }
        }
        else {
            plot(orders, zeta.val, xlab = "number of sites", 
                ylab = "zeta-diversity", pch = 20, ylim = c(0, 
                  zeta.val[1]))
            lines(orders, zeta.val)
        }
        plot(orders, zeta.val, log = "y", pch = 20, xlab = "Zeta order", 
            ylab = "Zeta-diversity", main = "Exponential regression")
        lines(orders, 10^predict.lm(zeta.exp, data.frame(orders)))
        plot(orders, zeta.val, log = "xy", pch = 20, xlab = "Zeta order", 
            ylab = "Zeta-diversity", main = "Power law regression")
        lines(orders, 10^predict.lm(zeta.pl, data.frame(orders)))
    }
    return(zeta)
}

# Apply Zeta.decline to the list of clusters (Warning, it takes long).

zetavalues <- lapply(clusters_ready,Zeta.decline.v2) # The function Zeta.decline.v2 is a modified version of
						     # the package function Zeta.decline() that allows to run it
						     # through lapply(). Run the code for Zeta.decline.v2 before 
						     # executing this line. 
save.image('zetavalues.RData')

