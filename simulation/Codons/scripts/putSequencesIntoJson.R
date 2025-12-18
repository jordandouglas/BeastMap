if (!require(rjson , quietly = T)) install.packages("rjson")
library(rjson)




# Load json
JSON = fromJSON(file="var.json")



# Protein sequences
seqFile = "simulated.xml"
xml = readLines(seqFile)
xml = gsub('"', "'", xml)

	

# Get simulated sequences
sequences = xml[grep("<sequence id", xml)]
sequences = sequences[grep("taxon[=]'", sequences)]
sequences = sequences[ (length(sequences)/2 + 1) : length(sequences)]
sequences = gsub("id[=].+spec", "spec", sequences)


str = paste(sequences, collapse="\n")
JSON[["seq-data"]] = str


# Read parameters
sim.df = read.table("simulated.log", header=T, sep = "\t")
sub.df = sim.df[1,]
for(x in colnames(sub.df)){
	val = sub.df[,x]
	JSON[[x]] = val
}	
		

JSON_str = as.character(rjson::toJSON(JSON, indent=1))
write(JSON_str, "var.seq.json")


