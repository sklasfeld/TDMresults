#get name of file
args <- commandArgs(trailingOnly = TRUE)
print(paste("File_Examined:",args[1]))

## examine the object via a connection, which will be opened as needed.
con <- gzfile(args[1])
print(readRDS(con))
#str(readRDS(con))
close(con)
