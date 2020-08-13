
####################################################################################################
rename_samples <- function(samples, session_root_dir){
    # ********************************************************************************************
    #
    # ********************************************************************************************

    # Load file that contains the sample renaming info.
    sample_renaming = load_sample_rename_file(
                                    input_file = file.path(session_root_dir, SAMPLE_RENAME_FILE),
                                    sample_names = samples)


    # Rename samples in cell and tissue segmentation files.
    # ****************************************************
    for(old_name in samples){
        new_name = sample_renaming[[old_name]]
        old_dir = file.path(session_root_dir, old_name)
        new_dir = file.path(session_root_dir, new_name)
        stopifnot(dir.exists(old_dir))
        stopifnot(!dir.exists(new_dir))
        for(x in new_dir) dir.create(x, recursive=FALSE)
        log_message(paste0('rename ', old_name, ' to ', paste(new_name, collapse=' / '),
                           ifelse(length(new_name)>1, ' [sample will be split]','')),
                    level=2)

        # Load cell and tissue segmentation files, change sample name and write to disk.
        for(file_extension in c(CELL_FILES_EXTENSION, TISSUE_FILES_EXTENSION)){
            old_file = file.path(old_dir, paste0(old_name, file_extension))
            new_file = file.path(new_dir, paste0(new_name, file_extension))
            tmp = read.table(old_file, sep='\t', as.is=TRUE, header=TRUE,
                             check.names=TRUE, strip.white=TRUE, colClasses='character')
            stopifnot('sample_name' %in% colnames(tmp))
            stopifnot('image_id' %in% colnames(tmp))

            # Case 1: the data does not need to be split.
            if(length(new_name) == 1){
                tmp[,'sample_name'] = new_name
                write.table(tmp, file=new_file, sep='\t', quote=FALSE, row.names=FALSE)

            # Case 2: the data must be split into 2.
            } else if(length(new_name) == 2){
                # Get imageID values that corresponds to the 2 sub-samples in the original sample.
                if(file_extension == CELL_FILES_EXTENSION){
                    stopifnot('cell_y_position' %in% colnames(tmp))
                    image_id_split = split_by_coordinate(coordinates = tmp[,'cell_y_position'],
                                                         image_ids   = tmp[,'image_id'])
                }
                #stopifnot(all(unique(tmp[,'image_id']) %in% unlist(image_id_split)))
                if(!all(unique(tmp[,'image_id']) %in% unlist(image_id_split))) raise_error(
                    msg = 'Tissue seg data has more image IDs than cell seg data.',
                    file = old_file, type = 'warning')

                # Temporary check for the split.
                if(FALSE){
                    jpeg(filename = file.path(session_root_dir, paste0(old_name, '_split.jpg')),
                         width = 700, height = 1000)
                        rand_sample = sample(1:nrow(tmp), ifelse(nrow(tmp) > 5000, 5000, nrow(tmp)))
                        plot(tmp$cell_x_position[rand_sample], tmp$cell_y_position[rand_sample],
                             col=ifelse(tmp[rand_sample, 'image_id'] %in% image_id_split[[1]],
                                        'darkgreen', 'darkred'))
                    dev.off()
                }

                # Split original sample based on imageID values.
                tmp_1 = tmp[tmp[,'image_id'] %in% image_id_split[[1]], ]
                tmp_2 = tmp[tmp[,'image_id'] %in% image_id_split[[2]], ]
                if(file_extension == CELL_FILES_EXTENSION) stopifnot(nrow(tmp_1)+nrow(tmp_2) == nrow(tmp))
                tmp_1[,'sample_name'] = new_name[1]
                tmp_2[,'sample_name'] = new_name[2]
                write.table(tmp_1, file=new_file[1], sep='\t', quote=FALSE, row.names=FALSE)
                write.table(tmp_2, file=new_file[2], sep='\t', quote=FALSE, row.names=FALSE)

            } else stop("ERROR: THIS CASE IS NOT IMPLEMENTED YET.")
            unlink(old_file)
        }
        unlink(old_dir, recursive=T)
    }



    # Rename samples in PARAMETERS_FILE.
    # *********************************
    file_name = file.path(session_root_dir, PARAMETERS_FILE)
    file_content = read_file_as_vector(file_name, ignore_comments=FALSE)

    # Find start and end position of 'samples:' paramter in the parameter file.
    start_position = grep('^samples:', file_content)
    stopifnot(length(start_position) == 1)
    parameter_position = grep(':', file_content)
    stopifnot(length(parameter_position) >= 2)
    if(start_position == tail(parameter_position, 1)){
        end_position = length(file_content)
    } else end_position = parameter_position[which(parameter_position == start_position) + 1] - 1

    # Replace old sample values with new values.
    if(start_position == 1){
        file_content = c('samples:',
                         as.character(unlist(sample_renaming)),
                         file_content[(end_position + 1):length(file_content)])

    } else if(end_position == length(file_content)){
        file_content = c(file_content[1:(start_position - 1)],
                         'samples:',
                         as.character(unlist(sample_renaming)))

    } else{
        file_content = c(file_content[1:(start_position - 1)],
                         'samples:',
                         as.character(unlist(sample_renaming)),
                         file_content[(end_position + 1):length(file_content)])
    }

    # Overwrite paramter input file.
    write(file_content, file=file_name, append=F, ncolumns=1)


    # Rename samples in THRESHOLDS_FILE.
    # *********************************
    file_name = file.path(session_root_dir, THRESHOLDS_FILE)
    if(file.exists(file_name)){
        # Load thrsholds file.
        tmp = read.table(file_name, sep='\t', as.is=TRUE, header=TRUE,
                         check.names=TRUE, strip.white=TRUE, colClasses='character')
        # Legacy format support: change sample name column header.
        if(colnames(tmp)[1] == 'sampleName') colnames(tmp)[1] = 'sample_name'
        stopifnot('sample_name' %in% colnames(tmp))
        sample_values = tmp[,'sample_name']

        # Replace sample names.
        for(old_name in samples){
            new_name = sample_renaming[[old_name]]
            if(length(new_name) == 1){
                sample_values[which(sample_values == old_name)] = new_name
            } else if(length(new_name) == 2){
                row_nb = which(sample_values == old_name)
                sample_values[row_nb] = new_name[1]

                # Duplicate thresholds for second new name of sample.
                sample_values = c(sample_values, rep(new_name[2], length(row_nb)))
                tmp = rbind(tmp, tmp[row_nb,])

            } else stop("ERROR: THIS CASE IS NOT IMPLEMENTED YET.")

        }
        tmp[,'sample_name'] = sample_values
        tmp = tmp[order(tmp[,'sample_name']), ]
        write.table(tmp, file=file_name, sep='\t', quote=FALSE, row.names=FALSE)
    }

}
####################################################################################################



####################################################################################################
load_sample_rename_file <- function(input_file, sample_names){
    # ********************************************************************************************
    #
    # ********************************************************************************************

    # Load file content by line. Lines starting with # are ignored.
    file_content = read_file_as_vector(input_file)
    if(length(file_content) < 2) raise_error(
        msg  = 'Sample renaming files must contain at least 2 lines: header + one sample.',
        file = input_file)

    # For each line, extract the correspondence between the old and new sample names.
    sample_renaming = list()
    for(line in file_content[-1]){
        line_elements = string_to_vector(line)
        if(!length(line_elements) %in% c(2,3)) raise_error(
            msg  = 'All rows of sample renaming files must have either 2 or 3 elements.',
            file = input_file)
        sample_name_old = line_elements[1]
        sample_name_new = line_elements[2:length(line_elements)]

        # Verify the values we are adding are not duplicated
        if(sample_name_old %in% names(sample_renaming)) raise_error(
            msg = paste0('Duplicated sample name: ', sample_name_old),
            file = input_file)
        for(x in sample_name_new){
            if(x %in% as.character(unlist(sample_renaming)) |
               any(duplicated(sample_name_new))) raise_error(
                msg = paste0('Duplicated new sample name: ', x),
                file = input_file)
        }

        # Replace any ' character in the new sample name with "_"
        sample_renaming[[sample_name_old]] = gsub("'", '-', sample_name_new)
    }

    # Verify no new name is the same as one of the old names.
    for(x in 1:length(sample_renaming)){
        sample_name_new = sample_renaming[[x]]
        for(new_name in sample_name_new){
            if(new_name %in% names(sample_renaming[-x])) raise_error(
                msg = paste0('New sample name is same as one of the old names: ', new_name),
                file = input_file)
        }
    }

    return(sample_renaming)
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

