## this file gets the expected md5sums from the MANIFEST file
## and calculated the md5sums for the other files
## and makes sure they are the same

library(data.table)



#' compare the expected md5sums (from manifest file )
#' against the computed md5sums
#' print to screen the result
#' 
#' @param manifest_file, full file path to the manifest file
#' @param md5sums_file, full file path to the computed md5sums
compare_md5sums <- function(manifest_file, md5sums_file){

    manifest_df = data.frame(fread(manifest_file))
    manifest_df = manifest_df[,c("File", "MD5")]
    manifest_df$File = gsub("Fastq/", "", manifest_df$File)

    computed_df = data.frame(fread(md5sums_file, header=F))
    colnames(computed_df) = c("MD5_computed", "File")

    # normalize dirextory structure
    computed_df$File = normalizePath(computed_df$File)
    dir_str = paste0(dirname(md5sums_file), "/")
    computed_df$File = gsub(dir_str, "", computed_df$File)

    compare_df = merge(manifest_df, computed_df, all = T)
    print("number mismatched:")
    print(sum(compare_df$MD5 != compare_df$MD5_computed))

    print("file mD5 mismatched:")
    print(compare_df[which(compare_df$MD5 != compare_df$MD5_computed),])

}

#
### read in arguements
args = commandArgs(trailingOnly=TRUE)
manifest_file = args[1]
md5sums_file = args[2]
compare_md5sums(manifest_file, md5sums_file)

