### Script to summarize DCC files
### 2020-10-12 Tyler Hether thether@nanostring.com



### #######################
### User-defined Parameters
### #######################

# Provide the dcc directory (absolute or relative path)
input_dcc_directory <- "DCC"
# Provide a prefix that will prepend the output files.
output_prefix <- "output"
# Number of processors to use.
n_processors <- 4

# End User-defined Parameters. Run the remaining code as source 


### ########
### Preamble
### ########

# Start fresh but retain user-defined parameters
rm(list=setdiff(ls(), c("input_dcc_directory", "output_prefix", "n_processors")))

# List of packages for session
.packages = c("plyr", "dplyr", "parallel", "data.table")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
null <- lapply(.packages, require, character.only=TRUE)

### ##########################
### Check input path and files
### ##########################

# Add a trailing / to directory name
if(!endsWith(input_dcc_directory, "/")){
  input_dcc_directory <- paste0(input_dcc_directory, "/")
}
# Confirm that directory exists
if(!dir.exists(input_dcc_directory)){
  stop(paste0("The input directory, ", input_dcc_directory, 
              " was not found. Please check path."))
}
# Get the dcc files in input_dcc_directory
dcc_files <- Sys.glob(paste0(input_dcc_directory, "*dcc"))
# Confirm that there is at least 1 dcc file in input_dcc_directory
if(length(dcc_files)<1){
  stop(paste0("The input directory, ", input_dcc_directory, 
              " was found but no files ending with .dcc were located.\n", 
              "Are these files compressed?"))
}

### #########
### Functions
### #########

process_dcc <- function(the_dcc){
  # ARGS:
  # the_dcc is a character with a valid dcc file path.
  # 
  # returns a list of length two providing the parsed data
  
  # Here are the lines
  the_lines <- readLines(the_dcc)
  
  # Here is the Sample_ID
  scan_attributes_header_line <- which(the_lines=="<Scan_Attributes>")
  if(length(scan_attributes_header_line)!=1){
    stop(paste0("The Scan_Attributes header was not found in ", the_dcc))
  }
  if(grepl("^ID,", the_lines[scan_attributes_header_line+1])){
    Sample_ID <- strsplit(the_lines[scan_attributes_header_line+1], split="ID,")[[1]][2]
  } else {
    stop(paste0("The Sample_ID line was not found in ", the_dcc))
  }
  
  # Raw, Trimmed, Stitched, and Aligned are bounded withing
  # the NGS_Processing_Attributes header
  # Here are the bounding indeces
  ngs_prcoessing_header_line <- which(the_lines=="<NGS_Processing_Attributes>")
  if(length(ngs_prcoessing_header_line)!=1){
    stop(paste0("The NGS_Processing_Attributes header was not found in ", the_dcc))
  }
  ngs_prcoessing_footer_line <- which(the_lines=="</NGS_Processing_Attributes>")
  if(length(ngs_prcoessing_footer_line)!=1){
    stop(paste0("The NGS_Processing_Attributes footer was not found in ", the_dcc))
  }
  # And parse the Raw, Trimmed, Stitched, and Aligned
  rtsa_names <- c("Raw", "Trimmed", "Stitched", "Aligned")
  rtsa <- lapply(rtsa_names, function(x){
    x_line <- which(grepl(paste0("^",eval(x), ","), 
                    the_lines[ngs_prcoessing_header_line:ngs_prcoessing_footer_line]))
    if(length(x_line)!=1){
      stop(paste0("The ", x, " counts were not found in file ", the_dcc))
    } else {
      x_out <- as.numeric(
        strsplit(the_lines[ngs_prcoessing_header_line:ngs_prcoessing_footer_line][x_line], 
               split=paste0(eval(x), ","))[[1]][2])
      return(x_out)
    }
  })
  rtsa_df <- as.data.frame(do.call(cbind, rtsa))
  colnames(rtsa_df) <- rtsa_names
  # Make summary_df by combining data (missing Uniques at this point)
  summary_df <- cbind(data.frame(Sample_ID=Sample_ID), rtsa_df)
  
  # Now pull out the RTS IDs and their associated counts.
  code_summary_header_line <- which(the_lines=="<Code_Summary>")
  if(length(code_summary_header_line)!=1){
    stop(paste0("The Code_Summary header was not found in ", the_dcc))
  }
  code_summary_footer_line <- which(the_lines=="</Code_Summary>")
  if(length(code_summary_footer_line)!=1){
    stop(paste0("The Code_Summary footer was not found in ", the_dcc))
  } 
  # If there are not entries, provide a warning
  # example:
  # <Code_Summary>
  # </Code_Summary>
  uniques <- 0
  if(code_summary_footer_line - code_summary_header_line < 2){
    warning(paste0("There are no entries in code summary for file ", the_dcc, 
                   " so all RTS values will be 0."))
    RTSs_df <- data.frame(Sample_ID=Sample_ID)
  } else {
    # Parse and transform
    RTSs <- the_lines[(code_summary_header_line+1):(code_summary_footer_line-1)]
    RTSs_long <- do.call(rbind, strsplit(RTSs, split=","))
    RTSs_wide <- matrix(as.numeric(RTSs_long[,2]), nrow=1, dimnames=list(NULL, RTSs_long[,1]))
    uniques <- apply(RTSs_wide, 1, sum)
    RTSs_df <- data.frame(Sample_ID=Sample_ID, as.data.frame(RTSs_wide))
  }
  
  # Add the uniqes to summary_df
  summary_df$Unique <- uniques
  
  # Return the two data.frame objects
  return(list(summary_df, RTSs_df))
  
}


### ##########
### Processing
### ##########

# Two steps: 1) Process each file 2) merge data

## Step 1 process each file
# Set up cluster
cl <- makeCluster(n_processors)
clusterExport(cl=cl, varlist=c("dcc_files", "process_dcc"), envir=environment())
# Execute the main function in parallel
all_processed <- parLapply(cl, dcc_files, process_dcc)
# Stop cluster
stopCluster(cl)

## Step 2 merge data
# The summary file.
summary_out <- do.call(rbind, lapply(all_processed, "[[", 1L))

# The dedup file. 
# Note that the RTS IDs need to be in the same order so
# we will use dplyr's bind_rows
dedup_out <- bind_rows(lapply(all_processed, "[[", 2L))

# Ensure there are no conversion issues with Sample_IDs.
if(any(is.na(dedup_out$Sample_ID))){
  stop("Sample_ID was converted to NA, implying an error in processing.")
}

# Let the user know that how many NA counts are being converted to zero and convert.
message(paste0("Converting ", length(which(is.na(dedup_out))), " NAs to zero."))
# dedup_out[is.na(dedup_out)] <- 0

# reorder the column names
cols_sorted <- c("Sample_ID", sort(colnames(dedup_out)[2:ncol(dedup_out)]))
dedup_out <- dedup_out %>% select(eval(cols_sorted))

# Why convert NA to zero? 
# When columns do not exist, bind_rows places an NA in the position.
# x <- data.frame("Sample_ID"="one", "RTS1"=1, "RTS2"=2)
# y <- data.frame("Sample_ID"="two", "RTS3"=3, "RTS2"=2.1, "RTS4"=4)
# z <- bind_rows(x,y)
# z
# z[is.na(z)] <- 0
# z

### ##########
### Write data
### ##########

message("Writing data")
data.table::fwrite(summary_out, file = paste0(output_prefix, "_summary.txt"), 
            sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
data.table::fwrite(x=dedup_out, file = paste0(output_prefix, "_dedup.txt"), 
                   sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, na=0)


message("End of script")
# End of script