if (!require(rjson , quietly = T)) install.packages("rjson")
library(rjson)


BURNIN = 0
NSIMS = 100


getTrees = function(fileName){

	file.in = readLines(fileName, warn=F)
	trees = file.in[grep("tree STATE_", file.in)]
	trees = gsub(".+[=] ", "", trees)
	
	# Translate
	trans = file.in[(grep("Translate", file.in)+1):(grep("^;$", file.in)-1)[1]]
	trans = gsub("\t", "", trans)
	trans = gsub("^ +", "", trans)
	trans = gsub(",", "", trans)
	indexes = sapply(strsplit(trans, " "), function(ele) ele[1])
	labels  = sapply(strsplit(trans, " "), function(ele) ele[2])
	
	
	for (i in 1:length(indexes)){
	
		trees = gsub(paste0("[(]", indexes[i], "[:]"), paste0("(", labels[i], ":"), trees)
		trees = gsub(paste0(",", indexes[i], "[:]"), paste0(",", labels[i], ":"), trees)
		trees = gsub(paste0("[(]", indexes[i], "[[]"), paste0("(", labels[i], "["), trees)
		trees = gsub(paste0(",", indexes[i], "[[]"), paste0(",", labels[i], "["), trees)
	
	}
	
	trees

}




# Read parameters
truth.df = read.table("truth/truth.log", header=T, sep = "\t")

# Read species tree
species.trees = getTrees("truth/truth.trees")



# Postprocessing burnin
burnin.start = floor(BURNIN*nrow(truth.df))
if (burnin.start <= 0) burnin.start = 1
include = seq(from = burnin.start, to = nrow(truth.df), length = NSIMS)
include = floor(include)



# Prepare folders
for (simnum in 1:NSIMS){

	rownum = include[simnum]
	
	
	#tmp
	#simnum = simnum - 1

	f = paste0("templates/rep", simnum)
	cat(paste(f, "\n"))
	dir.create(f, showWarnings=F) 
	
	
	
	# Create json
	JSON = list()
	
	
	# Prior variables
	sub.df = truth.df[rownum,]
	for(x in colnames(sub.df)){
		val = sub.df[,x]
		JSON[[x]] = val
	}	
	
	
	# Tree
	JSON[["tree"]] = gsub("&", "&amp;", species.trees[rownum])
	write(JSON[["tree"]], paste0(f, "/trueTree.newick"))

	
	# frequencies vector
	frequencies = as.numeric(truth.df[rownum,grep("frequencies", colnames(truth.df))])
	JSON[["frequencies"]] = paste0(frequencies, collapse = " ")


	JSON_str = as.character(rjson::toJSON(JSON, indent=1))
	write(JSON_str, paste0(f, "/var.json"))

}









