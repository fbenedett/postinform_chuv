
####################################################################################################
#' Rename samples in the input data and split them into two if needed.
#'
rename_samples <- function(sample_rename, root_dir){

    # Rename samples in cell and tissue segmentation files.
    # ****************************************************
    for(old_name in names(sample_rename)){
        new_name = sample_rename[[old_name]]
        old_dir = file.path(root_dir, old_name)
        new_dir = file.path(root_dir, new_name)
        stopifnot(dir.exists(old_dir))
        stopifnot(!dir.exists(new_dir))
        for(x in new_dir) dir.create(x, recursive=FALSE)
        log_message(paste0('rename ', old_name, ' to ', paste(new_name, collapse=' / '),
                           ifelse(length(new_name)>1, ' [sample will be split]', '')),
                    level=2)

        # Load cell and tissue segmentation files, change sample name and write to disk.
        for(file_extension in c(CELL_FILES_EXTENSION, TISSUE_FILES_EXTENSION)){
            old_file = file.path(old_dir, paste0(old_name, file_extension))
            new_file = file.path(new_dir, paste0(new_name, file_extension))
            tmp = read.table(old_file, sep='\t', as.is=TRUE, header=TRUE,
                             check.names=TRUE, strip.white=TRUE, colClasses='character')

            # Case 1: the data does not need to be split.
            if(length(new_name) == 1){
                tmp[,'sample_name'] = new_name
                write.table(tmp, file=new_file, sep='\t', quote=FALSE, row.names=FALSE)

            # Case 2: the data must be split into 2.
            } else if(length(new_name) == 2){
                # Get image ID values that corresponds to the 2 sub-samples in the original sample.
                if(file_extension == CELL_FILES_EXTENSION){
                    image_id_split = split_by_coordinate(coordinates = tmp[,'cell_y_position'],
                                                         image_ids   = tmp[,'image_id'])
                    #plot(tmp$cell_x_position[rand_sample], tmp$cell_y_position[rand_sample], col=as.factor(tmp$image_id[rand_sample]), pch=18)
                }
                if(!all(unique(tmp[,'image_id']) %in% unlist(image_id_split))) raise_error(
                    msg = c("Tissue seg data has more image IDs than cell seg data.",
                            "Image ID values no present in the cell seg data will be removed."),
                    file = old_file,
                    type = 'warning')

                # Temporary check for the split.
                if(FALSE){
                    jpeg(filename = file.path(root_dir, paste0(old_name, '_split.jpg')),
                         width = 700, height = 1000)
                        rand_sample = sample(1:nrow(tmp), ifelse(nrow(tmp) > 5000, 5000, nrow(tmp)))
                        plot(tmp$cell_x_position[rand_sample], tmp$cell_y_position[rand_sample],
                             col=ifelse(tmp[rand_sample, 'image_id'] %in% image_id_split[[1]],
                                        'darkgreen', 'darkred'))
                    dev.off()
                }

                # Split original sample based on image ID values.
                tmp_1 = tmp[tmp[,'image_id'] %in% image_id_split[[1]], ]
                tmp_2 = tmp[tmp[,'image_id'] %in% image_id_split[[2]], ]
                if(file_extension == CELL_FILES_EXTENSION) stopifnot(nrow(tmp_1)+nrow(tmp_2) == nrow(tmp))
                tmp_1[,'sample_name'] = new_name[1]
                tmp_2[,'sample_name'] = new_name[2]
                write.table(tmp_1, file=new_file[1], sep='\t', quote=FALSE, row.names=FALSE)
                write.table(tmp_2, file=new_file[2], sep='\t', quote=FALSE, row.names=FALSE)

            } else stop("cannot split sample in more than 2. THIS CASE IS NOT IMPLEMENTED YET.")
            unlink(old_file)
        }
        unlink(old_dir, recursive=T)
    }
}
####################################################################################################



####################################################################################################
#' Load a "sample rename" file from disk.
#'
#' @param input_file [string] path of the "sample rename" file to load.
#'
load_sample_rename_file <- function(input_file){

    # Load file content by line. Lines starting with # are ignored.
    file_content = read_file_as_vector(input_file, ignore_comments=TRUE,
                                       ignore_empty_line=TRUE, remove_quotes=TRUE)
    if(length(file_content) < 2) raise_error(
        msg  = 'Sample renaming files must contain at least 2 lines: header + one sample.',
        file = input_file)

    # For each line, extract the correspondence between the old and new sample names.
    sample_rename = list()
    for(line in file_content[-1]){
        line_elements = string_to_vector(line)
        if(!length(line_elements) %in% c(2,3)) raise_error(
            msg  = "All rows of 'sample_rename.txt' files must have either 2 or 3 elements.",
            file = input_file)
        old_name = line_elements[1]
        new_name = line_elements[2:length(line_elements)]

        # Add new old/new sample name pair to list.
        # Replace any ' character in the new sample name with "_".
        sample_rename[[old_name]] = gsub("'", '-', new_name)
    }

    return(sample_rename)
}
####################################################################################################


####################################################################################################
#'
#' @param sample_rename [list]
#'
check_sample_rename <- function(sample_rename, original_samples){

    # Verify that no old/new sample names are duplicated.
    check_for_duplicates(names(sample_rename),
                         error_message = 'Duplicated sample names found in file:',
                         error_file = SAMPLE_RENAME_FILE)
    check_for_duplicates(unlist(sample_rename, use.names=F),
                         error_message = 'Duplicated new sample names found in file:',
                         error_file = SAMPLE_RENAME_FILE)
    check_for_duplicates(c(names(sample_rename), unlist(sample_rename, use.names=F)),
                         error_message = 'New sample name has same name as old sample:',
                         error_file = SAMPLE_RENAME_FILE)

    # Verify that there is a 1:1 match between the sample_rename values and the original
    # sample names as passed in the input parameter file.
    check_for_difference(original_samples, names(sample_rename),
        error_message = "One or more samples have no matching value in the 'sample rename' file:",
        error_file = SAMPLE_RENAME_FILE)
    check_for_difference(names(sample_rename), original_samples,
        error_message = "One or more values in 'sample rename' file do not match any actual sample:",
        error_file = SAMPLE_RENAME_FILE)
}
####################################################################################################


####################################################################################################
check_for_duplicates <- function(values_to_check, error_message, error_file){
    duplicates = values_to_check[duplicated(values_to_check)]
    if(length(duplicates) > 0) raise_error(msg = error_message,
                                          items_to_list = duplicates,
                                          file = error_file)
}
####################################################################################################


####################################################################################################
check_for_difference <- function(set_1, set_2, error_message, error_file){
    if(length(setdiff(set_1, set_2)) > 0) raise_error(msg = error_message,
                                                      items_to_list = setdiff(set_1, set_2),
                                                      file = error_file)
}
####################################################################################################


####################################################################################################
split_by_coordinate <- function(coordinates, image_ids){
    # ********************************************************************************************
    #
    # ********************************************************************************************
    coordinates = as.numeric(coordinates)
    stopifnot(all(!is.na(coordinates)))

    # Find position of largest gap in coordinates: this is where the split must be done.
    sorted_coordinates = sort(coordinates)
    gap_size = diff(sorted_coordinates)
    gap_max_position = which(gap_size == max(gap_size))
    stopifnot(length(gap_max_position) == 1)
    coordinate_split_position = sorted_coordinates[gap_max_position] + max(gap_size)/2

    # Split image ID values at the position of the coordinate split.
    id_list_1 = unique(image_ids[which(coordinates <= coordinate_split_position)])
    id_list_2 = unique(image_ids[which(coordinates > coordinate_split_position)])
    stopifnot(length(id_list_1) > 0)
    stopifnot(length(id_list_2) > 0)
    if(length(intersect(id_list_1, id_list_2)) > 0) raise_error(
        msg = 'Split by coordinate failed: some image ID values overalp the split:',
        items_to_list = intersect(id_list_1, id_list_2))
    if(any(sort(unique(image_ids)) != sort(c(id_list_1, id_list_2)))) raise_error(
        msg = 'Split by coordinate failed: some image ID are missing after split.')

    # Return the two lists of image ID values.
    return( list('id_list_1'=id_list_1, 'id_list_2'=id_list_2) )
}
####################################################################################################

