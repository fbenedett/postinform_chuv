
####################################################################################################
load_session_parameters <- function(session_root_dir){
    # ********************************************************************************************
    # Load and verify parameter presence and their values.
    #
    # ********************************************************************************************

    # Read parameters file. Verify all required arguments are present.
    arg_values = read_parameters_file(file.path(session_root_dir, PARAMETERS_FILE))

    # Sample list.
    # ***********
    # Check that at least one sample was provided.
    if(! 'samples' %in% names(arg_values)) raise_error(msg="Parameter 'samples' is missing.",
                                                        file = PARAMETERS_FILE)
    samples = sort(unique(arg_values[['samples']]))
    if(length(samples) == 0) raise_error(msg="No 'samples' values in file.", file=PARAMETERS_FILE)

    # Check that sample names do not contain '[' or ']'. These are not allowed in sample names.
    if(any(grepl(pattern='\\[|\\]', x=samples))) raise_error(
        "Some sample names contain the non-authorized characters '[' or ']'.", file=PARAMETERS_FILE)

    # Verify that no sample name contains the string '_merge' or '-merge'. These are reserved for
    # 'merged' files.
    if(any(grepl(pattern='[_-]merge$', x=samples))) raise_error(
        msg = "Sample names ending in '_merge' or '-merge' are not allowed.", file=PARAMETERS_FILE)


    # Tissues list.
    # ************
    # Check that at least one tissue type was provided and that all tissues are part of the
    # authorized list.
    if(! 'tissues' %in% names(arg_values)) raise_error(msg="Parameter 'tissues' is missing.",
                                                        file = PARAMETERS_FILE)
    arg_values[['tissues']] = tolower(arg_values[['tissues']])
    tissues = unique(arg_values[['tissues']])
    if(length(tissues) == 0) raise_error(msg="No 'tissues' values in file.", file=PARAMETERS_FILE)

    nonauthorized_tissues = tissues[which(!tissues %in% AUTHORIZED_TISSUES)]
    if(length(nonauthorized_tissues) > 0) raise_error(
        msg = c(paste("One or more tissue types are not in the authorized list:",
                    paste(nonauthorized_tissues, collapse=', ')),
                "The list of authorized tissues is the following:"),
        file = PARAMETERS_FILE,
        items_to_list = AUTHORIZED_TISSUES)


    # Markers lists.
    # *************
    # Check that either phenotyped or scored markers were provided (or both).
    if(! any(c('phenotyped_markers','scored_markers') %in% names(arg_values))) raise_error(
        msg  = "Neither the 'phenotyped_markers' nor the 'scored_markers' parameter is present.",
        file = PARAMETERS_FILE)
    markers_phenotyped = unique(arg_values[['phenotyped_markers']])
    markers_scored = unique(arg_values[['scored_markers']])
    if(length(markers_phenotyped) + length(markers_scored) == 0) raise_error(
        msg  = "No 'phenotyped_markers' nor 'scored_markers' are present in the file.",
        file = PARAMETERS_FILE)
    if(any(markers_phenotyped %in% markers_scored)) raise_error(
        msg  = "Duplicated values between 'phenotyped_markers' and 'scored_markers'.",
        file = PARAMETERS_FILE)

    # Verify that marker values are in the authorized list.
    check_marker_is_authorized(markers_phenotyped, marker_type='phenotyped')
    check_marker_is_authorized(markers_scored, marker_type='scored')


    # Marker combinations.
    # *******************
    # Merge marker combinations with single markers used for phenotyping/scorring, and verify
    # the values are correct.
    individual_cell_types = c(paste0(c(markers_phenotyped, markers_scored), 'p'))
    if('all' %in% tolower(arg_values[['marker_combinations']])){
        marker_combinations = unlist(lapply(1:length(individual_cell_types),
                                       FUN=function(y) combn(individual_cell_types,
                                                             y,
                                                             FUN=function(x) paste(x, collapse='_'),
                                                             simplify=T)))
    } else marker_combinations = unique(arg_values[['marker_combinations']], individual_cell_types)
    marker_combinations = check_and_fix_marker_combination_values(
                                                            marker_combinations,
                                                            markers_phenotyped,
                                                            markers_scored,
                                                            file_values_are_from = PARAMETERS_FILE)

    # Threshold table.
    # ***************
    # If scored markers are listed in the input, a threshold file must be provided that contains
    # the threshold values for each scored marker.
    if(length(markers_scored) > 0){
        if(!file.exists(file.path(session_root_dir, THRESHOLDS_FILE))) raise_error(
            msg = paste0("No threshold file [", THRESHOLDS_FILE ,"] was found in the input data."))
        thresholds = load_thresholds_file(input_file = file.path(session_root_dir,THRESHOLDS_FILE),
                                          markers_scored = markers_scored,
                                          sample_names = samples,
                                          tissue_types = tissues,
                                          rewrite_input = TRUE)
    } else thresholds = character(0)


    # Cell compartment value.
    # **********************
    # If no cell compartment value is passed, assume default value.
    cell_compartment = ifelse('cell_compartment' %in% names(arg_values),
                              arg_values[['cell_compartment']], DEFAULT_CELL_COMPARTMENT)
    # Verify that the cell compartment value is in the authorized list.
    if(!cell_compartment %in% AUTHORIZED_COMPARTMENTS) raise_error(
        msg = c(paste0('[', cell_compartment, '], is not a valid cell compartment value'),
                'The list of authorized cell compartment values is:'),
        items_to_list = AUTHORIZED_COMPARTMENTS)


    # Phenotype confidence threshold.
    # ******************************
    # If no phenotype confidence threshold is passed, assume default value.
    phenotype_threshold = ifelse('phenotype_confidence_threshold' %in% names(arg_values),
                                 arg_values[['phenotype_confidence_threshold']],
                                 DEFAULT_PHENOTYPE_CONFIDENCE_THRESHOLD)
    # Verify that the phenotype confidence threshold value is in the range [0-100].
    if(!is.numeric(phenotype_threshold) |
       phenotype_threshold < 0 | phenotype_threshold > 100) raise_error(
           msg = paste0('The value passed to the [phenotype_confidence_threshold] argument',
                        'must be a number in the range [0-100].'))


    # Create list with all input arguments and their values.
    # *****************************************************
    input_parameters = list()
    input_parameters$session_root_dir    = session_root_dir
    input_parameters$samples             = samples
    input_parameters$tissues             = tissues
    input_parameters$markers_phenotyped  = markers_phenotyped
    input_parameters$markers_scored      = markers_scored
    input_parameters$marker_combinations = marker_combinations
    input_parameters$thresholds          = thresholds
    input_parameters$cell_compartment    = cell_compartment
    input_parameters$phenotype_threshold = phenotype_threshold
    return(input_parameters)
}
####################################################################################################



####################################################################################################
#' Read a Post-inForm parameter file and return a list of all arguments and values it contains.
#' Each argument can have one or more values, provided as a vector.
read_parameters_file <- function(input_file){

    # Load file content. Lines starting with a '#' are ignored.
    file_content = read_file_as_vector(input_file)

    # Parse file content.
    # ******************
    # A line containing a ':' character indicates that a new input parameter is starting on this
    # line. the value before the ':' is the name of the parameter, while all values comming after
    # the ':' are the values for the parameter and are stored in a vector.
    parameter_list = list()
    parameter_name = NULL
    for(line in file_content){
        # Check whether the line contains a parameter name.
        if(grepl(pattern='^[a-zA-Z]+[ _a-zA-Z]*:[^:]*$', x=line)){
            # If this a parameter has was already read, save it to the parameter list before
            # moving on to the next parameter.
            if(!is.null(parameter_name)){
                parameter_list[[parameter_name]] = string_to_vector(parameter_values)
            }

            # Set parameter name and values.
            tmp = trimws(unlist(strsplit(line, split=':')))
            stopifnot(length(tmp) == 1 | length(tmp) == 2)
            parameter_name = gsub(' ', '_', tmp[1])
            parameter_name = tolower(gsub('_[_]+', '_', parameter_name))
            parameter_values = ifelse(length(tmp) == 2, tmp[2], '')

        } else{
            # Verify that a parameter name has already been defined, and that the parameter value
            # line does not contain any ':' character.
            if(is.null(parameter_name)) raise_error(
                msg="The first non-comment line of parameter files must contain a parameter value.",
                file=input_file)
            if(grepl(pattern=':', x=line)) raise_error(
                msg="':' characters are only allowed as parameter delimitors.", file=input_file)

            # Append line content to parameter values.
            parameter_values = trimws(paste(parameter_values, line, sep=" "))
        }

    }
    # Add last parameter and its value to the list.
    parameter_list[[parameter_name]] = string_to_vector(parameter_values)
    return(parameter_list)
}
####################################################################################################



####################################################################################################
read_file_as_vector <- function(input_file, ignore_comments=TRUE, ignore_empty_line=TRUE){
    # ********************************************************************************************
    # Read a text file from disk and return its content as vector of strings, where each element
    # corresponds to a line in the file.
    # In addition, white spaces are trimmed, and lines starting with a # character are ignored.
    # ********************************************************************************************
    stopifnot(file.exists(input_file))
    file_connection = file(input_file, open='r', encoding=guess_file_encoding(input_file))
    file_content = readLines(con=file_connection)
    close(file_connection)
    file_content = trimws(file_content)

    if(ignore_empty_line) file_content = file_content[which(file_content != '')]
    if(ignore_comments)   file_content = file_content[which(!startsWith(file_content, '#'))]
    return(file_content)
}
####################################################################################################



####################################################################################################
string_to_vector <- function(input_string){
    # Split a string into a vector of strings by ',', ';', tabs or spaces.
    input_string = trimws(gsub(pattern=',|;|\t', replacement=' ', x=input_string))
    if(input_string=='') return(character(0))
    return( unlist(strsplit(input_string, split='[ ]+')) )
}
####################################################################################################



####################################################################################################
load_thresholds_file <- function(input_file,
                                 markers_scored,
                                 sample_names,
                                 tissue_types,
                                 rewrite_input=FALSE){
    # ********************************************************************************************
    # Input arguments:
    #   input_file    : string. Path + name of threshold file to load.
    #   markers_scored: string vector. List of markers that should be found in the input file.
    #   sample_names  : string vector. List of samples that should be found in the input file.
    #   rewrite_input : If TRUE, the original input file is replaced with the (possibly) curated
    #                   content produced by the function.
    # ********************************************************************************************

    # Load input file and standardize column names.
    # ********************************************
    # Load file content by line. Lines starting with # are ignored.
    file_content = read_file_as_vector(input_file)
    if(length(file_content) < 2) raise_error(
        msg  = "Thresholds files must contain at least 2 lines: header + one sample.",
        file = input_file)

    # Split lines into elements and convert to table.
    file_content = lapply(file_content, FUN=function(x) string_to_vector(x))
    if(length(unique(sapply(file_content, FUN=length))) > 1) raise_error(
        msg  = "Not all lines in the file have the same number of elements.",
        file = input_file)

    # Standardize column names.
    col_names = file_content[[1]]
    col_names[grep('^sample', col_names, ignore.case=T)] = 'sample_name'
    col_names[grep('^tissue', col_names, ignore.case=T)] = 'tissue_type'
    if(col_names[1] != 'sample_name') raise_error(
        msg  = "The first column of the threshold file must be 'sample_name'.",
        file = input_file)
    if('tissue_type' %in% col_names & col_names[2] != 'sample_name') raise_error(
        msg  = "If present, the 'tissue_type' column must be the second column of the file.",
        file = input_file)


    # Convert input to data frame.
    # ***************************
    file_content = file_content[-1]
    thresholds = data.frame(matrix(unlist(file_content), nrow=length(file_content), byrow=T),
                            stringsAsFactors=FALSE)
    colnames(thresholds) = col_names

    # Add a 'tissue_type' column if none was provided.
    if(!'tissue_type' %in% col_names){
        thresholds = data.frame('sample_name' = thresholds[,1],
                                'tissue_type' = rep(tissue_types, each=nrow(thresholds)),
                                thresholds[,2:ncol(thresholds)],
                                stringsAsFactors = FALSE)
        col_names = c('sample_name', 'tissue_type', col_names[-1])
        colnames(thresholds) = col_names
    }
    thresholds = thresholds[order(thresholds[,1]), ]

    # Set correct data type for each column.
    character_cols = col_names %in% c('sample_name','tissue_type')
    thresholds[,which(character_cols)]  = sapply(thresholds[,which(character_cols)], as.character)
    thresholds[,which(!character_cols)] = sapply(thresholds[,which(!character_cols)], as.numeric)
    if(any(is.na(thresholds))) raise_error(
        msg  = "Non-numeric threshold values detected in input file.",
        file = input_file)


    # Verify input data.
    # *****************
    # Check for duplicated columns and rows.
    if(any(duplicated(col_names))) raise_error(msg = 'Input file has duplicated column names.',
                                               file = threshold_file)
    if(any(duplicated(thresholds, MARGIN=1))) raise_error(msg = 'Input file has duplicated rows.',
                                                          file = threshold_file)

    # Check that all expected markers and samples are present in the data frame.
    for(x in c('markers', 'samples')){
        if(x=='markers') missing_values = markers_scored[which(!markers_scored %in% col_names)]
        if(x=='samples') missing_values = sample_names[which(!sample_names %in% thresholds[,1])]
        if(length(missing_values) > 0) raise_error(
            msg = paste0("One or more ", x, " are missing from the threshold file:"),
            file = input_file,
            items_to_list = missing_values)
    }

    # Verify and correct tissue type values.
    thresholds[, 2] = check_and_fix_tissue_type(tissue_type_values = thresholds[, 2],
                                                file_values_are_from = input_file)

    # Verify all samples vs. tissue type combinations exist.
    for(sample in sample_names){
    for(tissue in tissue_types){
        row_nb = which(thresholds[,'sample_name'] == sample & thresholds[,'tissue_type'] == tissue)
        if(length(row_nb) == 0) raise_error(
            msg = paste0("No threshold values found for sample [", sample,
                         "] and tissue type [",tissue,"]."),
            file = input_file)
        if(length(row_nb) > 1) raise_error(
            msg = paste0("More than one threshold value found for sample [", sample,
                         "] and tissue type [",tissue,"]."),
            file = input_file)
    }
    }


    # Re-write file on disk, if asked to do so.
    # ****************************************
    if(rewrite_input) write.table(thresholds, file=input_file, quote=F, sep='\t', row.names=F)
    return(thresholds)

}
####################################################################################################



####################################################################################################
scan_dir_for_seg_data_files <- function(input_directory){
    # ********************************************************************************************
    # Scan the input_directory for *_cell_seg_data and *_tissue_seg_data_summary files. Verify
    # the files go in pairs, and return a list with two elements: 'cell_files' and 'tissue_files'.
    #
    # ********************************************************************************************
    stopifnot(dir.exists(input_directory))

    # List of all *_cell_seg_data.txt and *_tissue_seg_data_summary.txt files in input directory.
    cell_files = list.files(path=input_directory, pattern=paste0('.*', CELL_FILES_EXTENSION, '$'),
                            all.files=FALSE, full.names=FALSE, recursive=FALSE, ignore.case=FALSE)
    tissue_files = list.files(path=input_directory, pattern=paste0('.*',TISSUE_FILES_EXTENSION,'$'),
                              all.files=FALSE, full.names=FALSE, recursive=FALSE, ignore.case=FALSE)

    # If no files were found, return NULL.
    if(length(cell_files) + length(tissue_files) == 0) return(NULL)

    # Verify that the cell and tissue segmentation files go by pairs.
    prefix_cell   = gsub(pattern='_cell_seg_data.txt', replacement='', x=cell_files)
    prefix_tissue = gsub(pattern='_tissue_seg_data_summary.txt', replacement='', x=tissue_files)
    prefix_cell   = sort(prefix_cell)
    prefix_tissue = sort(prefix_tissue)
    difference_list = c(setdiff(prefix_cell,prefix_tissue), setdiff(prefix_cell,prefix_tissue))
    if(length(difference_list) > 0) raise_error(
        msg  = "The following cell and/or tissue segmentation files are not in pairs:",
        file = input_directory,
        items_to_list = difference_list)

    # The function returns the prefix of all cell/tissue segmentation files.
    return(prefix_cell)

}
####################################################################################################
