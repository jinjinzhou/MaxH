## File
file <- "estHsq.txt"
 
## Create connection
con <- file(description=file, open="r")

linn <-readLines(con)
n.lines <- length(linn)
name.vec <- rep(NA,n.lines)
gvar.vec <- rep(NA,n.lines)
for (i in 1:n.lines){  
   tmp <- unlist(strsplit(linn[i], " ", fixed = TRUE))
   name.vec[i] = tmp[1]
   gvar.vec[i] = tmp[3]
}
close(con)

file <- "../datasets/CG10kNHWRes4Subtyping.txt"
con <- file(description=file, open="r")
line <- readLines(con,n=1)
names <- unlist(strsplit(line, " ", fixed = TRUE))

n.pheno <- length(names)-2
gvar.matrix <- matrix(0,n.pheno,n.pheno)

k=1
for (i in seq(1,n.pheno)){
	gvar.matrix[i,i] = as.numeric(gvar.vec[k])
	k=k+1
	for (j in seq(i+1,n.pheno)){
		if(j>n.pheno) break
		gvar.matrix[i,j] = as.numeric(gvar.vec[k])
		k=k+1
	}
}
close(con)
