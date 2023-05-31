
source("http://bioconductor.org/biocLite.R")
library(Biostrings)

setup_sample_sheet <- function(sheet)
{
	
	
	#sheet$Plate.No. = gsub(" ", "_", sheet$Plate.No.)
	sheet = as.data.frame(apply(sheet,2,function(x)gsub('\\s+', '',x)))
	sheet$Plate.No. = gsub("plateC", "plate_C", sheet$Plate.No.)
	sheet$Plate.No. = gsub("plate1", "plate_1", sheet$Plate.No.)
	sheet["i5.index.reverse.complement"] = reverse_comp(sheet$i5.index)
	sheet["Barcode"] = paste(sheet$i7.index, sheet$i5.index.reverse.complement, sep="_")
	sheet["well.type"] = get_well_type(sheet$Cell) 
	sheet["Full.Sample.Name"] = paste(sheet$cell.type, sheet$Sample.Name, sheet$Plate.No., sheet$well.type, sep="_")
	sheet$Full.Sample.Name = trim(sheet$Full.Sample.Name)
	return (sheet)
}

trim <- Vectorize(function(x)
{
	if(substr(x, nchar(x), nchar(x)) == "_")
	{
		return (substr(x, 1, nchar(x)-1))
	}else
	{
		return (as.character(x))
	}
})

get_well_type <- Vectorize(function(x)
{
	if(toupper(x) == "SINGLE")
	{
		return ("")
	}else
	{
		return (as.character(x))
	}
})

reverse_comp <- Vectorize(function(x)
{
	as.character(reverseComplement(DNAString(x)))
})
