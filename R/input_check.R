####################################################################################################
#' Scans the input directory and searches for the following files:
#'  - sub-directories containing "_cell_seg_data.txt" and "_tissue_seg_data_summary.txt" files.
#'  - PARAMETERS_FILE
#'  - SAMPLE_RENAME_FILE
#'
#' @param input_dir [string] Path of input directory to scan.
#' @param output_dir [string] Path of output directory.
#'
inputdir_check <- function(input_dir, output_dir){

    # Verify the specified input and output directories exist.
    if(!dir.exists(input_dir)) raise_error(paste('Input directory could not be found:', input_dir))
    if(!dir.exists(output_dir)) raise_error(paste('Output directory could not be found:', output_dir))

    # Verify that at least one sub-directory is present and that it contains at least one pair of
    # cell/tissue segmentation files. We do not check file content at this point.
    subdirs = list.dirs(input_dir, recursive=FALSE)
    if(length(subdirs) == 0 |
       all(sapply(subdirs,
                  FUN=function(x) length(scan_dir_for_seg_data_files(input_directory=x)) == 0))
       ) raise_error(
            msg = paste0('At least one subdirectory with *_cell_seg_data.txt and ',
                         '*_tissue_seg_data_summary.txt files must be present in the input data.'),
            file = input_dir)

    # Backward compatibility: if no 'parameters.txt' file is present, we attempt to generate one
    # from legacy "marker combination" and "threshold" files.
    file_list = list.files(path=input_dir, all.files=FALSE, full.names=FALSE, recursive=FALSE)
    file_list = file_list[!file.info(file.path(input_dir, file_list))$isdir]
    if(!PARAMETERS_FILE %in% file_list){
        if(any(grepl('^marker_combinations.*\\.txt$',
                     x=file_list))) generate_parameter_file_from_legacy_input(input_dir)
    }

    # Verify that input parameters and sample rename files are present, and if needed rename them.
    rename_file_by_pattern(file_name=PARAMETERS_FILE, pattern='param', dir_name=input_dir,
                           out_dir=output_dir, raise_error_if_absent=TRUE)
    rename_file_by_pattern(file_name=THRESHOLDS_FILE, pattern='threshold', dir_name=input_dir,
                           out_dir=output_dir, raise_error_if_absent=FALSE)
    rename_file_by_pattern(file_name=SAMPLE_RENAME_FILE, pattern='rename', dir_name=input_dir,
                           out_dir=output_dir, raise_error_if_absent=FALSE)
    return(invisible(NULL))
}
####################################################################################################



####################################################################################################
#' Search for a single file whose name matches @pattern in directory @dir_name. If a file is found,
#' but its name is not exactly @out_dir/@file_name, rename/copy it to @out_dir/@file_name.
#'
#' @param file_name [string]
#' @param raise_error_if_absent [logical] If TRUE, and error is raised if not exactly 1 file
#'     matching @pattern is found in @dir_name.
#'
rename_file_by_pattern <- function(file_name, pattern, dir_name,
                                   out_dir = dir_name,
                                   raise_error_if_absent = FALSE){

    # Verify input and output directories exist.
    if(!dir.exists(dir_name)) raise_error("Unable to find input directory.", file=dir_name)
    if(!dir.exists(out_dir)) raise_error("Unable to find output directory.", file=out_dir)

    # List all files in input directory.
    file_list = list.files(path=dir_name, all.files=FALSE, full.names=FALSE, recursive=F)
    file_list = file_list[!file.info(file.path(dir_name, file_list))$isdir]

    # If the target file already exists, copy it to the output directory, if needed.
    if(file_name %in% file_list){
        if(out_dir != dir_name) file.copy(from=file.path(dir_name, file_name), to=out_dir)
        return(invisible(NULL))
    }

    # Rename single file matching the search pattern.
    file_to_rename = grep(pattern, file_list, ignore.case=TRUE, value=TRUE)
    if(length(file_to_rename) == 1){
        if(out_dir == dir_name){
            file.rename(from=file.path(dir_name, file_to_rename), to=file.path(out_dir, file_name))
        } else{
            file.copy(from=file.path(dir_name, file_to_rename), to=file.path(out_dir, file_name))
        }
        return(invisible(NULL))
    }

    # If no or multiple matches are found, raise error.
    if(length(file_to_rename) > 1) raise_error(
        msg = paste0('Multiple potential [', file_name, '] files found in input directory:'),
        file = dir_name,
        items_to_list = file_to_rename)
    if(raise_error_if_absent & length(file_to_rename) == 0) raise_error(
        msg = paste0('No [', file_name, '] file found in input directory.'),
        file = dir_name)
    return(invisible(NULL))
}
####################################################################################################



####################################################################################################
standardize_and_split_cell_data <- function(input_file,
                                            samples,
                                            markers_phenotyped,
                                            markers_scored,
                                            cell_compartment,
                                            phenotype_confidence_threshold,
                                            delete_input_file = FALSE){
    # ********************************************************************************************
    # Differences between inForm versions:
    # version 2.2
    #  - *_tissue_seg_data_summary.txt files contain a column named "Region Area (pixels)".
    #
    # version 2.4
    #  - a new column named "Annotation ID" is added to *_cell_seg_data.txt.
    #  - *_tissue_seg_data_summary.txt files contain column named "Region Area (square microns)"
    #  - in addition the "Sample ID" column of the *_cell_seg_data.txt files no longer contains
    #    the image ID values. These are now present in the "Annotation ID" column.
    #
    #
    # Input arguments:
    #   input_file:
    #   markers_phenotyped: vector of strings. List of input markers that are phenotyped.
    #   markers_scored
    #   cell_compartment
    #   phenotype_confidence_threshold
    #   delete_input_file: if TRUE, the input_file is deleted after it was split by samples.
    # ********************************************************************************************

    # Load input table. Verify it is not empty and standardize the column names.
    input_table = read.table(input_file, sep='\t', as.is=T, header=T,
                             colClasses='character', check.names=T, strip.white=T)
    if(nrow(input_table) == 0) raise_error('Input file has zero rows.', file = input_file)
    colnames(input_table) = standardize_column_names(colnames(input_table), input_file)


    # Verify all required columns are present in the table. If the 'phenotype' and/or 'confidence'
    # columns are missing, dummy values are created for them.
    if(!'phenotype' %in% colnames(input_table)) input_table[,'phenotype'] = 'MISSING'
    if(!'confidence' %in% colnames(input_table)) input_table[,'confidence'] = 100
    verify_columns_are_present(required_columns = c('sample_name', 'tissue_category', 'phenotype',
                                                    'confidence', 'cell_id', 'cell_x_position',
                                                    'cell_y_position'),
                               actual_columns   = colnames(input_table),
                               file_name        = input_file)


    # Extract sample names and filter to keep only rows matching samples names present in the
    # input parameters.
    sample_names = extract_sample_name(input_table[,'sample_name'], input_file=input_file)
    input_table = input_table[sample_names %in% samples,]
    sample_names = sample_names[sample_names %in% samples]
    if(nrow(input_table) == 0) raise_error(
        msg = 'No matching values found in [sample_name] column for any of the input samples',
        file = input_file,
        items_to_list = samples)


    # Extract image ID values. If the "annotation_id" column is present (inForm 2.4), extract
    # image ID values from it. Otherwise extract from the "sample_name" column (inForm 2.2).
    colnb = tail(which(colnames(input_table) %in% c('sample_name', 'annotation_id')), n=1)
    image_ids = extract_imageid(input_table[,colnb], input_file=input_file)


    # Capitalize tissue types in the "tissue_category" column, if needed.
    tissue_type_values = check_and_fix_tissue_type(
        tissue_type_values = input_table[,'tissue_category'],
        file_values_are_from = input_file)


    # Verify and fix 'Phenotype' and 'Confidence' values.
    confidence_values = check_and_fix_confidence_values(
        confidence_values    = input_table[,'confidence'],
        file_values_are_from = input_file)
    phenotype_values = check_and_fix_phenotype_values(
        phenotype_values     = input_table[,'phenotype'],
        markers_phenotyped   = markers_phenotyped,
        file_values_are_from = input_file)


    # Re-class 'phenotype_values' for rows where confidence < phenotype_confidence_threshold
    # to the value of 'MISSING'.
    phenotype_values[which(confidence_values < phenotype_confidence_threshold)] = 'MISSING'


    # Verify marker intensity values. Verify all scored markers have matching '_mean' column.
    marker_column_list = which(endsWith(colnames(input_table), '_mean'))
    for(i in marker_column_list) input_table[,i] = check_and_fix_marker_intensity_values(
                                                                intensity_values = input_table[,i],
                                                                file_values_are_from = input_file)
    if(length(markers_scored) > 0){
        missing_markers = which(! paste0(markers_scored, '_mean') %in% colnames(input_table))
        if(length(missing_markers) > 0) raise_error(
            msg  = 'Missing scored marker columns for the following markers:',
            file = input_file,
            items_to_list = markers_scored[missing_markers])
    }


    # Split data by sample and save it to disk
    # ****************************************
    # Test whether file contains individual markers.
    prefix = paste0('individualmarker_', paste(markers_in_file_name(input_file), collapse=''), '_')
    if(prefix == 'individualmarker__') prefix = ''

    # Split data by samples.
    for(sample in sort(unique(sample_names))){
        rows_to_keep = which(sample_names==sample)
        out_dir = file.path(dirname(dirname(input_file)), sample)
        if(!dir.exists(out_dir)) dir.create(out_dir)
        out_file = file.path(out_dir, paste0(prefix, sample, CELL_FILES_EXTENSION))
        out_file_exists = file.exists(out_file)
        write.table(data.frame('sample_name'     = sample_names[rows_to_keep],
                               'image_id'        = image_ids[rows_to_keep],
                               'cell_x_position' = input_table[rows_to_keep, 'cell_x_position'],
                               'cell_y_position' = input_table[rows_to_keep, 'cell_y_position'],
                               'cell_id'         = input_table[rows_to_keep, 'cell_id'],
                               'tissue_category' = tissue_type_values[rows_to_keep],
                               'phenotype'       = phenotype_values[rows_to_keep],
                               input_table[rows_to_keep, marker_column_list],
                               stringsAsFactors  = FALSE),
                    file      = out_file,
                    append    = out_file_exists,
                    col.names = !out_file_exists,
                    sep       = '\t',
                    quote     = FALSE,
                    row.names = FALSE)
    }
    # Delete the original input file from disk.
    if(delete_input_file) unlink(input_file)

}
####################################################################################################



####################################################################################################
standardize_and_split_tissue_data <- function(input_file,
                                              samples,
                                              delete_input_file = FALSE){
    # ********************************************************************************************
    #
    # Input arguments:
    #  -> input_file:  string. Path + name of tissue data file to process.
    #  -> delete_input_file: logical. If TRUE, the original input_file is deleted at the end of
    #                        the function, after it was split by samples.
    #
    # Return value:
    # The function does not return anything. It writes the split files directly to disk.
    # ********************************************************************************************

    # Standardize and verify column names of input file.
    # *************************************************
    input_table = read.table(input_file, sep='\t', as.is=T, header=T,
                             colClasses='character', check.names=T, strip.white=T)
    if(nrow(input_table) == 0) raise_error('Input file has zero rows.', input_file = input_file)
    colnames(input_table) = standardize_column_names(colnames(input_table))
    verify_columns_are_present(required_columns = c('sample_name', 'tissue_category',
                                                    'region_area_percent'),
                               actual_columns   = colnames(input_table),
                               file_name        = input_file)


    # Extract sample name and "image ID" values.
    # *****************************************
    # Extract sample names and filter to keep only rows matching samples names present in the
    # input parameters.
    sample_names = extract_sample_name(input_table[,'sample_name'], input_file=input_file)
    input_table = input_table[sample_names %in% samples,]
    sample_names = sample_names[sample_names %in% samples]
    if(nrow(input_table) == 0) raise_error('Input file has zero rows for samples:',
                                           paste(samples, collapse=', '), file=input_file)


    # Extract image ID values. If the "annotation_id" column is present (inForm 2.4), extract
    # image ID values from it. Otherwise extract from the "sample_name" column (inForm 2.2).
    colnb = tail(which(colnames(input_table) %in% c('sample_name', 'annotation_id')), n=1)
    image_ids = extract_imageid(input_table[,colnb], input_file=input_file)


    # Capitalize tissue types in the "tissue category" column, if needed.
    # ******************************************************************
    tissue_type_values = check_and_fix_tissue_type(
        tissue_type_values = input_table[,'tissue_category'],
        file_values_are_from = input_file)


    # Extract tissue region surface values.
    # ************************************
    # Tissue surface values are recorded differently between inForm versions:
    #  -> inForm 2.2: surface is given in pixels and we thus have to convert it to square microns.
    #  -> inForm 2.4: surface is given directly in square microns.
    #
    # Case of inForm 2.4 files.
    if('region_area_square_microns' %in% colnames(input_table)){
        region_area_surface = as.numeric(input_table[, 'region_area_square_microns'])

        # Case of inForm 2.2 files.
        # Convert surfaces to square micrometers: 1 pixel is 0.246 square micrometers.
    } else if('region_area_pixels' %in% colnames(input_table)){
        region_area_surface = round(as.numeric(input_table[, 'region_area_pixels']) * 0.246, 1)

        # Non-recognized format. Raise error.
    } else{
        raise_error(msg = 'Missing [region_area_square_microns] or [region_area_pixels] column.',
                    file = input_file)
    }
    region_area_surface[which(is.na(region_area_surface))] = 0


    # Verify tissue region percentage values.
    # **************************************
    # Remove trailing '%' character, replace ',' by '.' (in case ',' was used as number fractional
    # part separator) and replace empty values with 0. Then convert values to number.
    region_area_percent = sub('%', '', input_table[, 'region_area_percent'])
    region_area_percent = sub(',', '.', region_area_percent)
    region_area_percent[which(region_area_percent=='')] = 0
    region_area_percent = as.numeric(region_area_percent)
    if(any(is.na(region_area_percent))){
        raise_error(msg = 'non-numeric or missing value in [Region Area (percent)] column.',
                    file = input_file)
    }


    # Split data by sample and save it to disk.
    # ****************************************
    # Test whether file contains individual markers.
    prefix = paste0('individualmarker_', paste(markers_in_file_name(input_file), collapse=''), '_')
    if(prefix == 'individualmarker__') prefix = ''

    # Split data by sample.
    for(sample in sort(unique(sample_names))){
        rows_to_keep = which(sample_names==sample)
        out_dir = file.path(dirname(dirname(input_file)), sample)
        if(!dir.exists(out_dir)) dir.create(out_dir)
        out_file = file.path(out_dir, paste0(prefix, sample, TISSUE_FILES_EXTENSION))
        out_file_exists = file.exists(out_file)
        write.table(sort_by_imageid_samplename_and_tissuecategory(data.frame(
                        'sample_name'         = sample_names[rows_to_keep],
                        'image_id'            = image_ids[rows_to_keep],
                        'tissue_category'     = tissue_type_values[rows_to_keep],
                        'region_area_surface' = region_area_surface[rows_to_keep],
                        'region_area_percent' = region_area_percent[rows_to_keep],
                        stringsAsFactors      = FALSE)),
                    file=out_file,
                    append=out_file_exists,
                    col.names=!out_file_exists,
                    sep='\t',
                    quote=FALSE,
                    row.names=FALSE)
    }
    # Delete the original input file from disk.
    if(delete_input_file) unlink(input_file)
}
####################################################################################################



####################################################################################################
#' Standardize column names of input files.
#'
#' @param column_names [string vector] Names of columns to standardize.
#' @param input_file [string] Path and name of file from which the columns were taken. Only used to
#'     display an error message.
#' @return Standardized column names.
standardize_column_names = function(column_names, input_file){

    # Replace any '.' in column names by an '_'. The '.' are generally introduced in column names
    # by R as a replacement of a non-authorized character such as a blank space or a bracket.
    # For readability, multiple '.' are replaced by a single '_'.
    column_names = gsub(pattern='\\.+', replacement='_', x=column_names)
    column_names = gsub(pattern='_+$', replacement='', x=column_names)

    # Convert standard column and cell compartment names to lower case.
    standard_names = c('sample_name', 'tissue_category', 'phenotype', 'cell_id',
                       'cell_x_position', 'cell_y_position', 'annotation_id', 'confidence',
                       'region_area_square_microns', 'region_area_percent', 'region_area_pixels')
    for(x in c(standard_names, AUTHORIZED_COMPARTMENTS)){
        column_names = sub(pattern=paste0('^', x), replacement=x, x=column_names, ignore.case=T)
    }

    # Convert '_Mean_' to '_mean_' (in "Mean_Normalized_Counts").
    column_names = sub(pattern='_mean_', replacement='_mean_', x=column_names, ignore.case=T)

    # Rename columns that contain mean marker intensity values.
    col_start_regexp = paste0('^(', paste(AUTHORIZED_COMPARTMENTS, collapse = '|'), ')_')
    for(i in grep(paste0(col_start_regexp, '.*_mean_.*'), x=column_names)){
        marker_name = sub(col_start_regexp, '', column_names[i])
        marker_name = sub('_.*$', '', marker_name)

        # If the marker is present in the AUTHORIZED_MARKERS list, correct its capitalization if
        # needed.
        x = which(toupper(marker_name) == toupper(AUTHORIZED_MARKERS))
        stopifnot(length(x) <= 1)
        marker_name = ifelse(length(x) == 0, marker_name, AUTHORIZED_MARKERS[x])

        # Rename column.
        column_names[i] = paste0(marker_name, '_mean')
    }

    # Verify there are no duplicate columns.
    duplicated_columns = which(duplicated(column_names))
    if(length(duplicated_columns) > 0) raise_error(
        msg  = 'Duplicated column names found in input file:',
        file = input_file,
        items_to_list = column_names[duplicated_columns])

    return(column_names)
}
####################################################################################################



####################################################################################################
check_and_fix_tissue_type = function(tissue_type_values, file_values_are_from){
    # ********************************************************************************************
    # Convert tissue type values to lower case and verify they are part of the authorized tissue
    # list.
    # tissue_type_values:   string vector. List of values to check.
    # file_values_are_from: string. Path + name of files from where the tissue type values
    #                       were taken. This is only used to display an error message.
    # ********************************************************************************************

    # Convert values to lower case. Replace missing tissue values with the default NOTISSUE value.
    tissue_type_values = tolower(trimws(tissue_type_values))
    tissue_type_values[tissue_type_values %in% NOTISSUES_SYNONYMS] = NOTISSUE

    # Verify that all values in tissue_type_values are part of the authorized list.
    offending_values = unique(tissue_type_values[!tissue_type_values %in% AUTHORIZED_TISSUES])
    if(length(offending_values) > 0) raise_error(
        msg = c(paste0('Invalid tissue type(s) found in input file: ',
                       paste0(offending_values, collapse = ', '), '.'),
                'The list of authorized tissue values is the following:'),
        items_to_list = AUTHORIZED_TISSUES,
        file = file_values_are_from)
    return(tissue_type_values)
}
####################################################################################################



####################################################################################################
#' Verifies that phenotype values are correct, and if needed, corrects them.
#'
#' @param phenotype_values
#' @param markers_phenotyped  : List of phenotyped markers for the current session.
#' @param file_values_are_from [string] Name of file from where phenotype values were taken. Only
#'     used to display error message.
check_and_fix_phenotype_values <- function(phenotype_values,
                                           markers_phenotyped,
                                           file_values_are_from = ''){

    # Replace all empty phenotype values by 'MISSING'. This is the placeholder for missing
    # phenotype values.
    phenotype_values = trimws(as.character(phenotype_values))
    row_ids = which(phenotype_values == '')
    if(length(row_ids) > 0){
        phenotype_values[row_ids] = 'MISSING'

        # Display warning indicating how many values where modified.
        missing_fraction = round(length(row_ids) / length(phenotype_values) * 100, 2)
        if(missing_fraction > 1.0) raise_error(
            msg=c('Empty rows found in [Phenotype] column of input file: ',
                paste0('Number of offending values: ', length(row_ids), ' [',missing_fraction,'%]'),
                'The values of these rows has been set to [MISSING].'),
            file = file_values_are_from,
            type = 'warning')
    }

    # Replace

    # Substitute '-' with '_' in Phenotype values. This is for the case where a '-' was used as
    # separator value instead of a '_'.
    phenotype_values = gsub(pattern='-', replacement='_', phenotype_values)

    # Substitute any trailing '+' with 'p' in the phenotype values.  E.g. 'CD8+' becomes 'CD8p',
    # or 'CD8+_CD11+' becomes 'CD8p_CD11p'.
    phenotype_values = gsub(pattern='\\+', replacement='p', phenotype_values)

    # Substitute any duplicated '_' with a single '_'. Remove leading or trailing '_'.
    phenotype_values = gsub(pattern='_[_]+', replacement='_', phenotype_values)
    phenotype_values = gsub(pattern='^_|_$', replacement='', phenotype_values)

    # To save time, we exit the function here if all phenotype values are correct. If this
    # is not the case, then the function carries-on and will try to correct those values.
    if(length(get_bad_markers(phenotype_values,
                              authorized_markers = markers_phenotyped,
                              ignored = IGNORED_PHENOTYPES)) == 0) return(phenotype_values)


    # Replace accepted NO_PHENOTYPE synonym values with NO_PHENOTYPE.
    phenotype_values[which(toupper(phenotype_values) %in%
                           toupper(c(NO_PHENOTYPE, NO_PHENOTYPE_SYNONYMS)))] = NO_PHENOTYPE

    # Replace accepted stroma and tumor synonyms with 'DAPIp' and 'CKp' respectively.
    for(x in 1:2) phenotype_values = gsub(
        pattern=paste0(rep(switch(x, AUTHORIZED_STROMA_VALUES, AUTHORIZED_TUMOR_VALUES), each=2),
                       c('', 'p'), collapse='|'),
        replacement=switch(x, 'DAPIp', 'CKp'), x=phenotype_values, ignore.case=TRUE)


    # Correct capitalization and add a 'p' suffix (for 'positive') to any marker missing it.
    bad_values = get_bad_markers(phenotype_values,
                                 authorized_markers=AUTHORIZED_MARKERS, ignored=IGNORED_PHENOTYPES)
    authorized_values = paste0(AUTHORIZED_MARKERS,'p')
    for(bad_value in bad_values){
        correct_value = authorized_values[which(tolower(authorized_values) == tolower(bad_value))]
        if(length(correct_value) == 0) correct_value = authorized_values[which(
                                    tolower(authorized_values) == paste0(tolower(bad_value),'p'))]
        if(length(correct_value) != 1) next
        for(x in 1:4){
            regexp_pattern = switch(x, paste0('^',bad_value,'$'),
                                       paste0('^',bad_value,'_'),
                                       paste0('_',bad_value,'_'),
                                       paste0('_',bad_value,'$'))
            replacement = switch(x, correct_value,
                                    paste0(correct_value,'_'),
                                    paste0('_',correct_value,'_'),
                                    paste0('_',correct_value))
            bad_ids = grep(regexp_pattern, phenotype_values, ignore.case = TRUE)
            if(length(bad_ids) > 0) phenotype_values[bad_ids] = gsub(pattern = regexp_pattern,
                                                                     replacement = replacement,
                                                                     x = phenotype_values[bad_ids],
                                                                     ignore.case=TRUE)
        }
    }
    stopifnot(sum(grepl('__', phenotype_values)) == 0)


    # At this point all values should be in the AUTHORIZED_MARKERS list.
    bad_values = get_bad_markers(phenotype_values,
                                 authorized_markers=AUTHORIZED_MARKERS, ignored=IGNORED_PHENOTYPES)
    if(length(bad_values) > 0) raise_error(
        msg = c(paste('One or more unauthorized markers were found in the [Phenotype] column:',
                      paste(bad_values, collapse=', ')),
                'The list of authorized markers is the following:'),
        file = file_values_are_from,
        items_to_list = AUTHORIZED_MARKERS)

    # If there are values not part of the phenotyped marker list, display a warning.
    bad_values = get_bad_markers(phenotype_values,
                                 authorized_markers=markers_phenotyped, ignored=IGNORED_PHENOTYPES)
    if(length(bad_values) > 0) raise_error(
        msg = 'Some marker values in the [Phenotype] column are not part of the input parameters:',
        file = file_values_are_from,
        items_to_list = bad_values, type = 'warning')

    return(phenotype_values)
}
####################################################################################################



####################################################################################################
get_bad_markers <- function(phenotype_values, authorized_markers, ignored){
    # ********************************************************************************************
    #
    # ********************************************************************************************
    marker_values = unique(unlist(strsplit(unique(phenotype_values), split='_')))
    authorized_list = c(paste0(authorized_markers,'p'), ignored)
    return(marker_values[which(! marker_values %in% authorized_list)])
}
####################################################################################################



####################################################################################################
check_and_fix_confidence_values <- function(confidence_values, file_values_are_from=''){
    # ********************************************************************************************
    # Takes a vector of confidence values, checks that they are correct, and if needed - and
    # possible - corrects them.
    #
    # Input arguments:
    #  -> confidence_values   : confidence values to check/fix.
    #  -> file_values_are_from: file name from which the confidence values are taken from.
    #                           This is only used to display in the error message.
    # ********************************************************************************************

    # Remove any trailing '%' character or whitespaces, and replace empty values '' with 0.
    confidence_values = gsub('%', '', trimws(confidence_values))
    confidence_values[which(confidence_values == '' | confidence_values == 'NaN')] = 0
    confidence_values = as.numeric(confidence_values)

    # Make sure that all values are now numeric.
    if(any(is.na(confidence_values))){
        offending_values = unique(confidence_values[which(is.na(confidence_values))])
        raise_error(msg = 'Non-numeric value detected in [Confidence] column:',
                    file = file_values_are_from,
                    items_to_list = offending_values)
    }
    return(confidence_values)
}
####################################################################################################


####################################################################################################
check_and_fix_marker_intensity_values <- function(intensity_values,
                                                  file_values_are_from = ''){
    # *******************************************************************************************
    # Takes a vector of marker intensity values, checks that they are correct, and if needed
    # corrects them.
    #
    # -> intensity_values    : marker intensity values to check/fix.
    # -> file_values_are_from: file name from which the marker intensity values are taken from.
    #                          This is only used to display in the error message.
    # *******************************************************************************************

    # Replace empty values with 0.
    intensity_values[which(grepl('^[[:space:]]*$', intensity_values))] = 0

    # Check whether values are in US format, i.e. they use a ',' instead of '.' as separator
    # between the integer and fractional part of the number.
    if(any(grepl(',', intensity_values))) intensity_values = sub(',','.',intensity_values)

    # Check for non-numeric values.
    tmp = unique(intensity_values)
    bad_ids = which(is.na(as.numeric(tmp)))
    if(length(bad_ids) > 0) raise_error(
        msg = 'Non-numeric values detected in one of the marker intensity columns:',
        file = file_values_are_from,
        items_to_list = tmp[bad_ids])

    # Convert all values to numeric.
    return(as.numeric(intensity_values))
}
####################################################################################################



####################################################################################################
check_and_fix_marker_combination_values <- function(marker_combinations,
                                                    markers_phenotyped,
                                                    markers_scored,
                                                    file_values_are_from = ''){
    # ********************************************************************************************
    # Verify marker combination (cell type) values.
    #   marker_list: list of accepted markers for the current session, i.e. list of phenotyped
    #                and scored markers for the current session.

    # ********************************************************************************************

    # Substitute '+' by 'p' if needed. Both '+' and 'p' can be used to indicate positive markers.
    # Remove duplicates.
    marker_combinations = gsub(pattern='\\+', replacement='p', x=marker_combinations)
    marker_combinations = unique(marker_combinations)

    # Check that all individual markers that compose a cell type (e.g. 'CD4p_CD8p' is decomposed in
    # 'CD4p' and 'CD8p') end in 'p' or 'm' (for positive/negative).
    individual_markers = unique(unlist(strsplit(marker_combinations, split='_')))
    if(any(! grepl(pattern="p$|m$", x=individual_markers) )) raise_error(
        'One or more cell types are missing the [p] or [m] suffix.', file=file_values_are_from)

    # Cell types that are composed of a single marker (e.g. "CD8p", "CD11p") must always end in
    # "p" (positive). Negative ("m") single markers make no sense.
    single_markers = marker_combinations[!grepl(pattern='_', x=marker_combinations)]
    if(any(! grepl(pattern="p$", x=single_markers))) raise_error(
        'one or more single cell types do not end in [p].', file=file_values_are_from)

    # Check that all individual markers belong to the authorized marker list.
    individual_markers = sub('p$|m$', '', individual_markers)
    check_marker_is_authorized(individual_markers, marker_type='combination')

    # Check that all individual markers are either part of the phenotyped or scored marker list.
    bad_values = individual_markers[!individual_markers %in% c(markers_phenotyped, markers_scored)]
    if(length(bad_values) > 0) raise_error(
        msg='The following markers present in combinations are not lised as phenotyped or scored.',
        items_to_list = bad_values,
        file = file_values_are_from)

    return(marker_combinations)
}
####################################################################################################



####################################################################################################
check_marker_is_authorized <- function(marker_list, marker_type){
    # ********************************************************************************************
    # Verify that all makers passed in the marker_list are part of the AUTHORIZED_MARKERS list,
    # and raise an error if not.
    #
    # Input arguments:
    #   marker_list: string vector. List of markers to check.
    #   marker_type: string. Type of marker, e.g. 'phenotyped', 'scored', 'combination'. This is
    #                used in the error message that is displayed.
    # ********************************************************************************************
    bad_markers = marker_list[!marker_list %in% AUTHORIZED_MARKERS]
    if(length(bad_markers) > 0) raise_error(
        msg = c(paste0('Non-authorized marker detected in [', marker_type, '] markers input: ',
                     paste(bad_markers, collapse= ', ')),
                'The list of authorized markers is the following:'),
        file = PARAMETERS_FILE,
        items_to_list = AUTHORIZED_MARKERS)
}
####################################################################################################



####################################################################################################
#' Test whether a file contains "individual marker" values. This is simply done by looking at the
#' file name. The convention is that individual marker files will contain at least one marker name
#' in their file name.
#' The function returns the list of markers in the file name, or character(0) if none is found.
markers_in_file_name <- function(file_name){
    return(sort(as.character(names(unlist(
        sapply(AUTHORIZED_MARKERS, FUN=function(x) grep(x, basename(file_name), ignore.case=T)))))))
}
####################################################################################################



####################################################################################################
#' Performs the following checks on the input data given a list of samples:
#'  * Input samples are present in the raw data.
#'  * Samples present in the raw data are also present in input parameter file [warning].
#'  * Each sample sub-directory contains exactly one cell and one tissue segmentation file.
#'
#' @param samples [str vector] sample names passed in input by the user.
#' @param root_dir [str] Root directory for the current session.
#' @return [NULL]
check_sample_data <- function(samples, root_dir){

    # Retrieve samples present in data based on sub-directory names present in session's root_dir.
    samples_from_data = sort(basename(list.dirs(root_dir, full.names=TRUE, recursive=F)))

    # Verify all samples passed in input by the user are present in the raw data.
    if(length(setdiff(samples, samples_from_data)) > 0) raise_error(
        msg = 'One or more samples listed in the input parameters are missing from the raw data:',
        items_to_list = setdiff(samples, samples_from_data))

    # Verify all samples present in the raw data are also present in the input parameters. If not,
    # raise a warning.
    if(length(setdiff(samples_from_data, samples)) > 0) raise_error(
        msg = 'One or more samples present in the data are not listed in input parameters:',
        items_to_list = setdiff(samples_from_data, samples),
        type = 'warning')

    # Verify each sample sub-directory contains exactly one cell and one tissue segmentation file.
    error_msg = NULL
    for(sample_name in samples){
        sample_dir = file.path(root_dir, sample_name)
        for(file_type in c(CELL_FILES_EXTENSION, TISSUE_FILES_EXTENSION)){
            files = list.files(path=sample_dir, pattern=paste0('.*', file_type, '$'),
                               all.files=F, full.names=F, recursive=F, ignore.case=F)

            if(length(files) == 0){
                error_msg = c(error_msg, paste0('Sample [', sample_name, ']: no [',
                                                file_type, '] file found.'))
            } else if(length(files) > 1){
                error_msg = c(error_msg, paste0('Sample [', sample_name,']: more than one [',
                                                file_type, '] file in subdirectory.'))
            }
        }
    }
    if(length(error_msg) > 0) raise_error(msg = "One or more errors detected in input data:",
                                          items_to_list = error_msg)
}
####################################################################################################



####################################################################################################
extract_sample_name <- function(input_vector, input_file='file not specified'){
    # ********************************************************************************************
    # Extract "sample name" values from input_vector strings. Depending on the version of inForm
    # that was used to generate the data, the format of the sample name is different:
    #  -> inForm 2.2: <sample name>[image id].im3, e.g.: OE37_V0[66030,14035].im3
    #  -> inForm 2.4: <sanple name>_Scan[number].qptiff, e.g. Bern_TMA1_Scan1.qptiff
    #
    # Input parameters:
    #  -> input_vector: vector of strings from which to extract sample names.
    #  -> input_file:   name of file from which the input strings were loaded from. This is only
    #                   needed to display an error message.
    # ********************************************************************************************

    # Define regexp for each supported version of inForm.
    regexp_inform22 = '\\_[[0-9]+,[0-9]+\\](_M[0-9])?\\.im3$'
    regexp_inform24 = '_Scan[0-9]+\\.qptiff$'


    # Case 1. Data format is inForm 2.2.
    inform22_pattern = grepl(pattern=regexp_inform22, x=input_vector, ignore.case=F)
    if(all(inform22_pattern)){
        return(sub(pattern=regexp_inform22, replacement='', x=input_vector))
    }
    # Case 2. Data format is inForm 2.4.
    inform24_pattern = grepl(pattern=regexp_inform24, x=input_vector, ignore.case=F)
    if(all(inform24_pattern)){
        return(sub(pattern=regexp_inform24, replacement='', x=input_vector))
    }

    # Case 3. Unsupported data format.
    if(all(inform22_pattern | inform24_pattern)){
        raise_error(msg='Mix of inForm 2.2. and 2.4 formats detected in [Sample Name] column.',
                    file=input_file)
    }
    if(sum(inform22_pattern) > 0){
        error_examples = input_vector[which(!inform22_pattern)]
    } else if(sum(inform24_pattern) > 0){
        error_examples = input_vector[which(!inform24_pattern)]
    } else{
        error_examples = input_vector
    }
    if(length(error_examples) > 10) error_examples = error_examples[1:10]
    raise_error(msg = 'Unsupported format detected in [Sample Name] column of input file:',
                file = input_file,
                items_to_list= error_examples)
}
####################################################################################################



####################################################################################################
sort_by_imageid_samplename_and_tissuecategory <- function(input_table){
    # ********************************************************************************************
    # Sort the input_table by its 'image_id' column values. Image ID values are composed of two
    # numbers separated by an '_', and represent coordinates of the image_id within the total
    # image file.
    # One problem is that sometimes the coordinates have 4 digits, and sometimes they have 5, but
    # using a simple sort on strings, '47430_12886' and '47430_8512' will be sorted in this order
    # (alphanumeric), when in fact it should be in the opposite order: "47430_8512", "47430_12886".
    # To solve that problem, the image ID strings are decomposed in individual coordinates, and
    # a sort based on the first, with a subsort based on the second is performed. Since these sorts
    # are now numeric and not alphabetical, things get put in the correct order.
    # ********************************************************************************************

    # Verify input data.
    stopifnot(is.data.frame(input_table))
    stopifnot('image_id' %in% colnames(input_table))
    stopifnot('sample_name' %in% colnames(input_table))
    stopifnot('tissue_category' %in% colnames(input_table))
    stopifnot(all(input_table[,'tissue_category'] %in% AUTHORIZED_TISSUES))

    # Split image ID values into their 2 components. Each image ID string must split in exactly
    # 2 values.
    tmp_list = strsplit(input_table$image_id, '_')
    stopifnot(all(sapply(tmp_list, length) == 2))

    # Retrieve individual coordinates. Make sure they are numeric values.
    coord1 = as.numeric(sapply(tmp_list, FUN=function(x) x[[1]]))
    coord2 = as.numeric(sapply(tmp_list, FUN=function(x) x[[2]]))
    stopifnot(all(!is.na(coord1)))
    stopifnot(all(!is.na(coord2)))

    # Attribute order value to tissue categories based on their oder in AUTHORIZED_TISSUES.
    tissue_category_order = as.numeric(sapply(input_table[,'tissue_category'],
                                              FUN=function(x) which(x == AUTHORIZED_TISSUES)))

    # Return the input table sorted by image ID values.
    input_table = input_table[order(input_table$sample_name,
                                    coord1, coord2, tissue_category_order, decreasing=F), ]
    row.names(input_table) = 1:nrow(input_table)
    return(input_table)
}
####################################################################################################

