####################################################################################################
#'
#'
#' @param input_dir
#' @param cell_compartment
#' @param output_dir
reduce_file_size <- function(input_dir, cell_compartment, output_dir = NULL){

    # Define columns to keep for cell and tissue segmentation files respectively.
    cols_to_keep_cell = c("Sample Name", "Tissue Category", "Phenotype", "Cell ID",
                          "Cell X Position", "Cell Y Position", "Confidence", "Annotation ID",
                          paste0(cell_compartment, ".* Mean .*"))
    cols_to_keep_tissue = c("Sample Name", "Tissue Category", "Region Area .*", "Annotation ID")


    # Loop through all files present in the input directory and reduce their size.
    #TODO: parallelized loop with "foreach" and "doParallel" libraries.
    # core_nb = parallel::detectCores() - 2
    # if(core_nb <= 0) core_nb = 1
    # if(core_nb > 4) core_nb = 4
    # cl = parallel::makeCluster(core_nb)
    # doParallel::registerDoParallel(cl)
    # foreach(prefix = scan_dir_for_seg_data_files(input_dir)) %dopar% {
    #     ...
    # }
    # parallel::stopCluster(cl)
    for(prefix in scan_dir_for_seg_data_files(input_dir)){
        cell_seg_file   = paste0(prefix, CELL_FILES_EXTENSION)
        tissue_seg_file = paste0(prefix, TISSUE_FILES_EXTENSION)

        log_message(file.path(basename(input_dir), cell_seg_file), level=2)
        delete_cols_in_file(input_file = file.path(input_dir, cell_seg_file),
                            colnames_to_keep = cols_to_keep_cell,
                            output_file = file.path(output_dir, cell_seg_file),
                            chunk_size = 200000,
                            sep = '\t')

        log_message(file.path(basename(input_dir), tissue_seg_file), level=2)
        delete_cols_in_file(input_file = file.path(input_dir, tissue_seg_file),
                            colnames_to_keep = cols_to_keep_tissue,
                            output_file = file.path(output_dir, tissue_seg_file),
                            chunk_size = 200000,
                            sep = '\t')
    }
    return(invisible(NULL))
}
####################################################################################################



####################################################################################################
#' Reads the input_file table chunk_size lines at a time and keep only columns listed in the
#' colnames_to_keep string vector.
#'
#' @param input_file [string] File to process.
#' @param colnames_to_keep [string vector]  Names of columns from input_file to keep.
#'     Regexp notation is allowed.
#' @param output_file [string] File where to save the output data table with the reduced number of
#'     columns. By default the output file is set to input_file, meaning that the input file gets
#'     overwritten.
#' @param chunk_size [integer] integer. Number of lines of the input file to read at a time. Larger
#'     values result in faster processing but larger memory consumption.
#' @param sep [string] Delimiter (separator) used in the input table file.
delete_cols_in_file <- function(input_file,
                                colnames_to_keep,
                                output_file = input_file,
                                chunk_size = 100000,
                                sep = '\t'){

    # Read file header and identify which columns of the file should be kept.
    stopifnot(file.exists(input_file))
    file_connection = file(input_file, open='r', encoding=guess_file_encoding(input_file))
    header = unlist(strsplit(readLines(con=file_connection, n=1), split=sep))
    cols_to_keep = grep(pattern=paste0('^', paste(colnames_to_keep, collapse='$|^'), '$'),
                        x=header, ignore.case=TRUE)

    # If there are no columns to remove, exit the function immediately to save time.
    if(length(cols_to_keep) == length(header)){
        close(file_connection)
        if(output_file != input_file){
            if(! dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive=T)
            file.copy(from=input_file, to=output_file, overwrite=FALSE)
        }
        return(invisible(NULL))
    }

    # Load content of file chunk_size lines at a time and keep only the subset of needed columns.
    reduced_df = NULL
    while(length(fchunk <- readLines(con=file_connection, n=chunk_size, warn=FALSE)) > 0){
        reduced_df = rbind(reduced_df,
                           as.data.frame(matrix(unlist(strsplit(fchunk, split=sep)),
                                                nrow=length(fchunk), byrow=T)[, cols_to_keep]))
        #message("read lines:", nrow(reduced_df))
    }
    close(file_connection)

    # Write the reduced data frame to output_file.
    colnames(reduced_df) = header[cols_to_keep]
    if(! dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive=T)
    write.table(reduced_df, file=output_file, row.names=FALSE, quote=FALSE, sep=sep)
    return(invisible(NULL))
}
####################################################################################################



####################################################################################################
#' Delete all files in the session_root_dir input directory that are not needed for the  analysis.
#'
delete_unnecessary_files <- function(session_root_dir){

    # Verify the input directory exists.
    if(!file.exists(session_root_dir) || !file.info(session_root_dir)$isdir) raise_error(
        paste('Input sub-directory could not be found:', session_root_dir))

    # Delete unnecessary files from all sub-directories present in the session's root dir.
    for(subdir in list.dirs(session_root_dir, full.names=TRUE, recursive=FALSE)){

        # Scan sub-directory for input files.
        # **********************************
        # Get list of all files in the current sub-directory. Delete the sub-directory if empty.
        file_list = list.files(path=subdir, all.files=FALSE, full.names=FALSE, recursive=FALSE)
        if(length(file_list) == 0){
            unlink(subdir, recursive=TRUE)
            next
        }

        # Identify 'rejected', 'seg_data' and 'merged' files. 'seg_data' files correspond to
        # *_cell_seg_data.txt and *_tissue_seg_data_summary.txt files.
        is_rejected = grepl(pattern='_rejected_', x=file_list, ignore.case=F)
        is_seg_data = grepl(pattern='.*_tissue_seg_data_summary.txt$|.*_cell_seg_data.txt$',
                               x=file_list, ignore.case=F)
        is_merge    = grepl(pattern='merge[_-]', x=file_list, ignore.case=T)


        # Delete un-necessary files.
        # *************************
        # Delete all files that are not needed for the analysis.
        # This includes:
        #  -> any file that is not a *_cell_seg_data.txt or *_tissue_seg_data_summary.txt file.
        #  -> *_rejected_* files. In principle only present if 'merge' files are present.
        #  -> if 'merge' files are present, any file that is not a merge file.

        # If 'merge' files are present, all 'non-merge' files get added to the list of files to
        # delete.
        if(any(is_merge)){
            to_delete = is_rejected | ! is_seg_data | ! is_merge
        } else to_delete = is_rejected | ! is_seg_data

        # Delete selected files, or the entire sub-directory if all files are to be deleted.
        if(all(to_delete)){
            unlink(subdir, recursive=T)
        } else sapply(file.path(subdir, file_list[which(to_delete)]), unlink)
    }
}
####################################################################################################



####################################################################################################
move_input_files_to_single_directory <- function(session_root_dir){
    # ********************************************************************************************
    # Move all cell and tissue segmentation files to single directory named SEGMENTATION_DATA_DIR.
    # ********************************************************************************************
    # Create the new 'seg_data' directory if needed.
    unique_dir = file.path(session_root_dir, SEGMENTATION_DATA_DIR)
    if(!dir.exists(unique_dir)) dir.create(unique_dir)

    # Move files in each directory to the new directory and delete the now empty directory.
    for(subdir in list.dirs(session_root_dir, full.names=TRUE, recursive=FALSE)){
        if(subdir==unique_dir) next
        lapply(list.files(subdir, all.files=FALSE, full.names=TRUE, recursive=F, include.dirs=F),
               FUN=function(x) file.rename(x, file.path(unique_dir, basename(x))) )
        unlink(subdir, recursive=T)
    }
}
####################################################################################################

