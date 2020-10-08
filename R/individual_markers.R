
####################################################################################################
merge_individual_marker_files <- function(input_dir){
    # ********************************************************************************************
    # Merge individual marker files located in input_dir into a single file. The merged result is
    # written to a new file and the original files are deleted from disk.
    #
    # ********************************************************************************************
    sample_name = basename(input_dir)
    for(extension in c(CELL_FILES_EXTENSION, TISSUE_FILES_EXTENSION)){

        # Identify files to merge. Make sure at least 2 files were found.
        files_to_merge = list.files(path=input_dir, pattern=paste0(extension, '$'),
                                    all.files=F, full.names=T, recursive=F, ignore.case=F)
        if(length(files_to_merge) == 1) raise_error(
            msg = paste('Individual marker merge error. Only one file to merge for sample:',
                        sample_name),
            file = input_dir)
        if(any(!startsWith(basename(files_to_merge), 'individualmarker_'))) raise_error(
            msg = 'Individual marker merge error: mix of individual and non-individual files.',
            file = input_dir)

        # Merge *_cell_seg_data.txt or *_tissue_seg_data_summary.txt files.
        if(extension==CELL_FILES_EXTENSION)   merged_df = merge_cell_data_files(files_to_merge)
        if(extension==TISSUE_FILES_EXTENSION) merged_df = merge_tissue_data_files(files_to_merge)

        # Write merged file to disk.
        file_name = file.path(input_dir, paste0(sample_name, extension))
        if(file.exists(file_name)) raise_error(
            msg = 'Individual marker merge error: file already exists.', file = file_name)
        write.table(merged_df, file=file_name, quote=FALSE, row.names=FALSE, sep='\t')

        # Delete original individual marker files.
        unlink(files_to_merge)
    }
}
####################################################################################################



####################################################################################################
merge_cell_data_files <- function(files_to_merge){
    # ********************************************************************************************
    # Load the provided *_cell_seg_data.txt files and merge their content into a single data
    # frame.
    #
    # Input arguments:
    #  -> files_to_merge : list[string], path + name of *_cell_seg_data.txt files to be merged.
    #
    # ********************************************************************************************
    stopifnot(length(files_to_merge) > 1)

    # Define columns that will be used as key fields during the merge.
    key_fields = c('sample_name', 'image_id', 'cell_x_position', 'cell_y_position')

    # Loop through all files that match the prefix + suffix combination. Load and check the
    # content of each file, then add it to the merged data frame.
    merged_df = NULL
    for(f in files_to_merge){

        # Load data and remove duplicated rows.
        input_df = remove_duplicated_rows(
                        input_df = read.table(f, sep='\t', header=TRUE, as.is=TRUE,
                                              colClasses='character', strip.white=TRUE),
                        key_fields = key_fields,
                        file_name = f)

        # Verify all input files have the same columns. This must always be the case as we have
        # already standardized the files earlier.
        if(is.null(merged_df)){
            stopifnot(all(colnames(input_df)[1:length(key_fields)] == key_fields))
            col_names = colnames(input_df)
        }
        stopifnot(all(colnames(input_df) == col_names))


        # Merge data.
        # **********
        # Merge data frame for the current marker with the global dataframe 'merged_df'.
        if(is.null(merged_df)){
            merged_df = input_df
        } else{
            # Append a suffix to the non-key fields, so they have a unique name and get preserved
            # during the merge.
            colnames(input_df)[(length(key_fields)+1):ncol(input_df)] = paste0(
                colnames(input_df)[(length(key_fields)+1):ncol(input_df)],
                '_', which(files_to_merge==f))

            row_count = nrow(merged_df)
            merged_df = merge(merged_df, input_df, by=key_fields, all=FALSE, sort=FALSE)

            # Compute percentage of data loss due to merge operation. If it's more than
            # FILE_LOSS_PERCENTAGE_THRESHOLD percent, a warning is raised.
            row_loss_percentage = round((row_count - nrow(merged_df)) / row_count * 100, 2)
            if(row_loss_percentage > FILE_LOSS_PERCENTAGE_THRESHOLD) raise_error(
                msg = c(paste0('Large number of records lost while merging cell data files [',
                               row_loss_percentage, '%].'),
                     'This could indicate a problem in the input data and should be investigated.'),
                file = f,
                type = 'warning')
        }
        rm(input_df)
    }


    # Extract tissue category, phenotype, and marker intensity values.
    # ***************************************************************
    # Create separate data frames for the non-key fields.
    tissue_cat_df = merged_df[, grep('tissue_category', colnames(merged_df))]
    phenotype_df  = merged_df[, grep('phenotype', colnames(merged_df))]
    marker_int_df = merged_df[, grep('_mean$|_mean_', colnames(merged_df))]
    for(x in 1:ncol(marker_int_df)) marker_int_df[,x] = as.numeric(marker_int_df[,x])

    # Remove the extracted columns from the merged data frame.
    merged_df = merged_df[,which(colnames(merged_df) %in% c(key_fields, 'cell_id'))]


    # Combine data from all "Tissue Category" columns.
    # ***********************************************
    # In principle, the "Tissue Category" columns of all individual markers should contain the
    # same value. Here we verify that this is the case and then keep only one copy of them.
    if(any(tissue_cat_df != tissue_cat_df[,1])){
        diff_rows = sort(unique(unlist(lapply(2:ncol(tissue_cat_df),
                                FUN=function(x) which(tissue_cat_df[x] != tissue_cat_df[,1])))))
        for(x in diff_rows){
            value_frequency = sort(table(as.character(tissue_cat_df[x,])), decreasing=T)
            stopifnot(length(value_frequency) >= 2)

            # Case 1: one of the tissue values has a majority within the row. Majority ruling is
            # possible and the most frequent value is used.
            if(as.numeric(value_frequency[1]) > as.numeric(value_frequency[2])){
                tissue_cat_df[x,] = names(value_frequency)[1]
                if(SHOW_TISSUE_CATEGORY_MISMATCH_WARNING) raise_error(
                    msg = c(paste0('[tissue_category] values differ across files. ',
                                   'Values were reconciled based on majority ruling.'),
                            paste0('Offending row: ', x)),
                    file = files_to_merge[1],
                    type = 'warning')

            # Case 2: majority ruling is not possible.
            } else{
                raise_error(
                    msg = c('Could not merge individual marker files.',
                            'Reason: tissue_category values differ across files with no majority.',
                            paste0('Offending row: ', x),
                            paste0('Offending values: ', paste(tissue_cat_df[x,], collapse=' '))),
                    file = files_to_merge[1])
            }
        }
        stopifnot(all(tissue_cat_df == tissue_cat_df[,1]))
    }
    # Add tissue category values to the merged data frame.
    merged_df[,'tissue_category'] = tissue_cat_df[,1]


    # Combine data from all 'Phenotype' columns.
    # *****************************************
    # Each "Phenotype" column contains the phenotyping for one or more individual marker. Here we
    # want to merge all these columns into a single column that will contain all phenotyping
    # information. If more than one marker is found for a given row, they get combined with a '_'
    # separator. E.g. if 'CD8p' and 'PD1p' are found in same row, the combination becomes
    # 'CD8p_PD1p'.
    # "Phenotype" values that are part of IGNORED_PHENOTYPES are removed at this point if they are
    # combined with other markers, so for instance 'CD8p_DAPIp_PD1p' becomes 'CD8p_PD1p'. If
    # all individual files contain "DAPIp", then a single "DAPIp" value is kept, e.g.
    # "DAPIp_DAPIp_DAPIp" becomes "DAPIp".
    regexp = paste(sapply(IGNORED_PHENOTYPES, FUN=function(x) paste0(x,'_|_',x)), collapse = '|')
    merged_df[,'phenotype'] = apply(phenotype_df, MARGIN=1,
                                    FUN=function(x) gsub(pattern = regexp,
                                                         replacement = '',
                                                         x = paste(unique(x), collapse='_')))


    # Combine data from all marker intensity columns.
    # **********************************************
    # Verify that all marker intensity columns in the merge data frame have the same values and
    # keep only one copy of them. A difference in values <= MARKER_INTENSITY_THRESHOLD is accepted.
    marker_intensity_cols = grep('_mean$', col_names, value=TRUE)
    for(col_name in marker_intensity_cols){
        col_index = grep(col_name, colnames(marker_int_df))
        differences = abs(marker_int_df[,col_index] - marker_int_df[,col_index[1]])
        differing_values = which(differences > MARKER_INTENSITY_THRESHOLD)

        # If there are values > MARKER_INTENSITY_THRESHOLD, an warning is displayed to the user.
        if(length(differing_values) > 0){
            differing_values = as.vector(as.matrix(differences))[differing_values]
            raise_error(
                msg=c(paste0('Values for column [', col_name, '] differ accross individual files'),
                      paste0('to merge by more than ', MARKER_INTENSITY_THRESHOLD, ' at ',
                             length(differing_values), ' occurences [',
                             round(length(differing_values)/nrow(marker_int_df)*100, 3), '%].'),
                      'Values from the first file (alphabetically) will be used.'),
                file=dirname(files_to_merge[1]),
                type = 'warning')
        }
    }

    # Add marker intensity values to merged dataframe.
    merged_df = cbind(merged_df, marker_int_df[, 1:length(marker_intensity_cols)])
    return(merged_df)
}
####################################################################################################



####################################################################################################
merge_tissue_data_files <- function(files_to_merge){
    # ********************************************************************************************
    # "merges" the content of all provided *_tissue_seg_data_summary.txt files into a single
    # data frame. Most of the time, all *_tissue_seg_data_summary.txt files are exactly the same,
    # but it can happen that some are missing a number of rows (i.e. some image subsets are
    # missing because there were excluded based on their poor quality).
    # So the "merge" operation done in this functions consists keeping only rows from the input
    # *_tissue_seg_data_summary.txt files that are found in all files (i.e. do the intersection
    # of all files.).
    #
    # ********************************************************************************************
    stopifnot(length(files_to_merge) > 1)

    # Define columns that will be used as key fields to carry out the merge.
    key_fields = c('sample_name', 'image_id', 'tissue_category')
    non_key_fields = c('region_area_surface', 'region_area_percent')

    # Merge all files in list.
    merged_df = NULL
    for(f in files_to_merge){

        # Load data and remove duplicated rows. Verify that the column names are correct - at this
        # point the files are standardized so they must all have the same column names.
        input_df = remove_duplicated_rows(
                        input_df = read.table(f, sep='\t', header=TRUE, as.is=TRUE, strip.white=T,
                                              colClasses=c(rep('character',3), rep('numeric',2))),
                        key_fields = key_fields,
                        file_name = f)
        stopifnot(all(colnames(input_df) == c(key_fields, non_key_fields)))

        # Merge data frame for the current marker with the global data frame 'merged_df'.
        if(is.null(merged_df)){
            merged_df = input_df
        } else{
            # Append suffix to non-key fields, so they have a unique name and get preserved
            # during the merge.
            colnames(input_df)[4:5] = paste0(non_key_fields, '_', which(files_to_merge %in% f))
            row_count = nrow(merged_df)
            merged_df = merge(merged_df, input_df, by=key_fields, all=FALSE, sort=FALSE)

            # Compute percentage of data loss due to merge operation. If it's more than
            # FILE_LOSS_PERCENTAGE_THRESHOLD percent, a warning is raised.
            row_loss_percentage = round((row_count - nrow(merged_df)) / row_count * 100, 2)
            if(row_loss_percentage > FILE_LOSS_PERCENTAGE_THRESHOLD) raise_error(
                msg = c(paste0('Large number of records lost while merging tissue data files [',
                               row_loss_percentage, '%].'),
                        'This could indicate a problem in the input data and should be investigated.'),
                file = f,
                type = 'warning')
        }
    }

    # Merge tissue surface values (absolute value or percentage).
    # The merge is made by keeping the smallest surface value from the individual files. This is
    # because, when merging cells value, we keep the intersection of all individual files, which,
    # in terms of surface, corresponds (roughly) to the smallest surface value in the tissue
    #surface files.
    for(i in 1:nrow(merged_df)){
        # Identify column with the smallest surface value.
        col_index = grep('region_area_surface', colnames(merged_df))
        min_index = which(merged_df[i, col_index] == min(merged_df[i, col_index]))[1]
        # Replace tissue surface values with the minimum value across all files.
        merged_df[i, 'region_area_surface'] = merged_df[i, col_index[min_index]]
        merged_df[i, 'region_area_percent'] = merged_df[i, col_index[min_index] + 1]
    }

    # Remove duplicated columns.
    return(merged_df[,1:5])
}
####################################################################################################



####################################################################################################
remove_duplicated_rows <- function(input_df, key_fields, file_name){
    # ********************************************************************************************
    # Check whether there are duplicated rows for key_fields in the input table (input_df), and if
    # so, removes them and displays a warning.
    #
    # file_name: name of file from where input_df is loaded. Only used to display warning to user.
    # ********************************************************************************************

    # Remove perfectly duplicated rows without displaying any warning-
    if(any(duplicated(input_df))) input_df = unique(input_df)

    # Remove rows duplicated over the key_fields with a warning to the user.
    duplicated_rows = which(duplicated(input_df[,key_fields]))
    if(length(duplicated_rows) > 0){
        input_df = input_df[-duplicated_rows,]
        raise_error(msg  = 'The following duplicated rows were deleted from input file:',
                    file = file_name,
                    items_to_list = duplicated_rows,
                    type = 'warning')
    }
    return(input_df)
}
####################################################################################################

