
####################################################################################################
generate_parameter_file_from_legacy_input <- function(session_root_dir){

    # ********************************************************************************************
    # Generate a 'parameters.txt' parameter file from the legacy 'marker_combinations_*.txt' and
    # 'thresholds_*.txt' files.
    # ********************************************************************************************

    # Get file name for 'marker_combinations_*.txt'
    marker_combination_file = list.files(path=session_root_dir,
                                         pattern='^marker_combinations.*\\.txt$',
                                         all.files=FALSE, full.names=TRUE,
                                         recursive=FALSE, ignore.case=FALSE)
    if(length(marker_combination_file) == 0) raise_error('No [marker_combinations_*.txt] in input.')
    if(length(marker_combination_file) > 1) raise_error('Multiple [marker_combinations_*.txt].')
    marker_combination_file = marker_combination_file[1]

    # Get file name for 'thresholds_*.txt'
    threshold_file = list.files(path=session_root_dir, pattern='^thresholds.*\\.txt$',
                                all.files=FALSE, full.names=TRUE,
                                recursive=FALSE, ignore.case=FALSE)
    if(length(threshold_file) == 0) raise_error('No [thresholds_*.txt] file in input.')
    if(length(threshold_file) > 1) raise_error('Multiple [thresholds_*.txt] files in input.')
    threshold_file = threshold_file[1]


    # Fix potential problems in 'marker_combinations_*.txt' and 'thresholds_*.txt' files.
    # Verify that file does not contain any spaces as column delimiters. We do this by checking
    # that the file doesn't contain any "space-tab", "tab-space" or "tab-tab" combination.
    # Lines are not allowed to start with space or tab.
    #if(fix_column_separator_in_table(marker_combination_file) == FALSE) stop('error detected')
    #if(fix_column_separator_in_table(threshold_file) == FALSE) stop('error detected')


    # Load cell and tissue types from 'marker_combinations_*.txt' file.
    tmp = load_marker_combination_data(marker_combination_file)
    cell_types = tmp$cell_types
    tissue_types = tmp$tissue_types

    # Load sample names and thresholds from 'thresholds_*.txt' file.
    thresholds = load_threshold_file_legacy(threshold_file, tissue_types)
    samples = as.character(dimnames(thresholds)[[1]])
    markers = as.character(dimnames(thresholds)[[2]])

    # Find which markers are phenotyped (threshold value of -1) and which are scored (threshold
    # value >= 0).
    phenotyped_markers = NULL
    scored_markers = NULL
    for(marker in markers){
        if(all(thresholds[,marker,] == -1)){
            phenotyped_markers = c(phenotyped_markers, marker)
        } else scored_markers = c(scored_markers, marker)
    }

    # Rename thresholds file.
    file.rename(threshold_file, file.path(session_root_dir, THRESHOLDS_FILE))

    # Create new "parameters.txt" file.
    file_connection = file(file.path(session_root_dir, "parameters.txt"), open="w")
    writeLines(paste0("# Auto-generated parameter file for ", basename(session_root_dir), "."),
               con=file_connection)
    writeLines("samples:", con=file_connection)
    for(x in samples) writeLines(x, con=file_connection)
    writeLines(paste0("tissues: ", paste(tissue_types, collapse=', ')), con=file_connection)
    writeLines(paste0("phenotyped_markers: ", paste(phenotyped_markers, collapse=', ')),
               con=file_connection)
    writeLines(paste0("scored_markers: ", paste(scored_markers, collapse=', ')),
               con=file_connection)
    writeLines("marker_combinations:", con=file_connection)
    for(x in cell_types) writeLines(x, con=file_connection)
    close(file_connection)
}
####################################################################################################



####################################################################################################
load_marker_combination_data <- function(marker_combination_file){

    # ********************************************************************************************
    # Reads an input 'marker_combinations_*.txt' file and returns a list with two elements:
    #  - cell_types  : a vector giving all cell types found in the input file.
    #  - tissue_types: a vector giving all tissue types found in the input file.
    #
    # Example of the function's return value:
    #  $cell_types
    #  [1] "CD4p_CD8p"     "CD8p_FOXP3p"   "KERATINp_CD8p"
    #  $tissue_types
    #  [1] "Stroma" "Tumor"
    #
    # Note: The cell type list is returned in the same order as the cell types are listed in the
    # input 'marker_combinations_*.txt' file.
    #
    # The combination file contains the list of cell types and the list of tissue types that we are
    # interested in for the current autodesc session. A cell type is a combination of markers, each
    # of which can be either positive "p" or negative "m". E.g. "CD4p_CD8p" means that the cell
    # type is positive for CD4 and positive for CD8, "CD8p_FOXP3m" means that the cell type is
    # positive for CD8 and negative for FOXP3. The list of cell types is found in the 'type' column
    # of the combination table. Example of input 'marker_combinations_*.txt' file.
    #   type	            tissues
    #   CD3p	            Stroma_Tumor
    #   CD3p_CD4p	        Stroma_Tumor
    #   CD8p_CD3p  	    Stroma_Tumor
    #   CD3p_CD4p_FOXP3p	Stroma_Tumor
    #   CD3p_CD8p_FOXP3p	Stroma_Tumor
    #
    # Input arguments:
    #   marker_combination_file: string, the full path+name of the combination file to load.
    # ********************************************************************************************

    # Load marker combination table and perform basic checks.
    # ******************************************************
    #  -> at least one row is present
    #  -> no duplicated rows or columns.
    #  -> table must contain the columns 'type' and 'tissues'.
    marker_combination_table = read.table(marker_combination_file, header=TRUE, sep='\t',
                                          as.is=TRUE, colClasses='character',
                                          fileEncoding=guess_file_encoding(marker_combination_file),
                                          strip.white=TRUE)
    file_name = basename(marker_combination_file)
    if(nrow(marker_combination_table) == 0) raise_error(msg  = 'Input file has zero rows.',
                                                        file = marker_combination_file)
    if(any(duplicated(colnames(marker_combination_table)))) raise_error(
        'Duplicated column name(s) in input file.', file=marker_combination_file)
    for(x in c('type', 'tissues')) if(! x %in% colnames(marker_combination_table)) raise_error(
        paste0('[', x, '] column missing in file.'), file=marker_combination_file)


    # Load marker combination values (cell types).
    # *******************************************
    # Make a list of all "cell types" listed in the input file. Each cell type must be unique, so
    # we check for duplicates.
    if(any(duplicated(marker_combination_table$type))) raise_error(
        'Duplicated cell type detected in input file.', file=marker_combination_file)
    cell_types = marker_combination_table$type

    # Substitute '+' by 'p' if needed. Both '+' and 'p' can be used to indicate positive markers.
    cell_types = gsub(pattern='\\+', replacement='p', x=cell_types)

    # Check that all individual markers that compose a cell type (e.g. 'CD4p_CD8p' is decomposed in
    # 'CD4p' and 'CD8p') end in "p" or "m" (for positive/negative).
    tmp = unique(unlist(strsplit(cell_types, split='_')))
    if(length(tmp) == 0) raise_error('No marker found in file.', file=marker_combination_file)
    if(any(! grepl(pattern='p$|m$', x=tmp) )) raise_error(
        'One or more cell types are missing the [p] or [m] suffix.', file=marker_combination_file)

    # Cell types that are composed of a single marker (e.g. "CD8p", "CD11p") must always end in
    # "p" (positive). Negative ("m") single markers make no sense.
    tmp = cell_types[!grepl(pattern='_', x=cell_types)]
    if(any(! grepl(pattern='p$', x=tmp))) raise_error(
        'One or more single cell types do not end in [p]', file=marker_combination_file)


    # Load tissue type values.
    # ***********************
    # Make a list of all tissue types (e.g. Stroma, Tumor) that are found in the
    # 'marker_combinations_*.txt' file. Tissue names must follow a fixed nomenclature and can only
    # have values listed in the AUTHORIZED_tissue_types, which is a global variable in autodesc.
    tissue_types = sort(unlist(strsplit(unique(marker_combination_table$tissues), split='_')))
    if(length(tissue_types) == 0) raise_error('No tissues found in input file.',
                                             file=marker_combination_file)

    # Verify tissue category values, and if needed, capitalize them.
    tissue_types = check_and_fix_tissue_type(tissue_type_values = tissue_types,
                                             file_values_are_from = file_name)


    # Return a list with two elements: the cell_types and the tissue_types.
    return(list('cell_types'=cell_types, 'tissue_types'=tissue_types))
}
####################################################################################################



####################################################################################################
load_threshold_file_legacy <- function(threshold_file, tissue_types){

    # ********************************************************************************************
    # Reads the input threshold_file passed by the user and returns a table with the thresholds
    # associated to each marker, tissue type and sample.
    # The threshold table contains the thresholds to be used to reclassify the intensity values
    # of each marker. 'intensity value >= threshold' determines whether a given cell belongs to
    # a cell type or not.
    # The threhold table can have multiple rows, where each row corresponds to a sample and
    # tissue type combination. The first column of the table must contain the sample name, the
    # second column is optional and contains the tissue type. All the subsequent colunms are the
    # thresholds for each marker. The column names are the names of the markers. Note that if
    # the 'tissue_type' column is missing, then autodesc assumes that the threshold values are
    # the same for all tissue types. Example of input threshold file:
    #
    #      sample_name  tissue_type  CD3  CD4  CD8  FOXP3
    #     SAMPLE_00801      Stroma   -1   -1   -1    1.1
    #     SAMPLE_00801 	     Tumor   -1   -1   -1	 1.1
    #     SAMPLE_00810 	    Stroma   -1   -1   -1	 0.9
    #     SAMPLE_00810       Tumor   -1   -1   -1	 0.9
    #
    #
    # Input arguments:
    #  -> threshold_file: string. the full path + name of the threshold file to load.
    #  -> tissue_types   : vector of strings. The list of tissue types that are considered
    #                     in the current autodesc session.
    # Return value:
    #  -> The function returns a 3D array, as illustrated here:
    #        , , Stroma
    #                     CD3 CD4 CD8 FOXP3
    #        SAMPLE_00801  -1  -1  -1  1.10
    #        SAMPLE_00810  -1  -1  -1  0.90
    #        SAMPLE_00816  -1  -1  -1  0.32
    #        , , Tumor
    #                     CD3 CD4 CD8 FOXP3
    #        SAMPLE_00801  -1  -1  -1  1.10
    #        SAMPLE_00810  -1  -1  -1  0.90
    #        SAMPLE_00816  -1  -1  -1  0.32
    #
    # ********************************************************************************************

    # Load table and perform basic checks:
    # ***********************************
    #  -> table must have at least one row
    #  -> no duplicated rows or columns.
    #  -> The first column of the table must be 'sample_name'.
    input_table = read.table(threshold_file, header=T, sep='\t', as.is=T,
                             fileEncoding=guess_file_encoding(threshold_file),
                             strip.white=T, colClasses='character')
    if(nrow(input_table) == 0) raise_error('Input file has zero rows.', file=threshold_file)
    col_names = colnames(input_table)
    if(any(duplicated(col_names))) raise_error('Input file has duplicated column names.',
                                               file=threshold_file)
    if(any(duplicated(input_table, MARGIN=1))) raise_error('Input file has duplicated rows.',
                                                           file=threshold_file)
    if(col_names[1] != 'sampleName') raise_error('first column of file must be named [sampleName].',
                                                 file=threshold_file)


    # Get list of samples from the threshold table.
    # ********************************************
    # Check sample names do not contain '[' or ']'. These are not allowed in sample names.
    sample_list = unique(input_table$sampleName)
    if(any(grepl(pattern='\\[|\\]', x=sample_list))) raise_error(
        "Some sample names contain the non-authorized characters '[' or ']'.", file=threshold_file)
    # Verify that no sample name contains the string '_merge' or '-merge'.
    # These are reserved for actual 'merged' files.
    if(any(grepl(pattern='[_-]merge$', x=sample_list))) raise_error(
        "Sample names ending in '_merge' or '-merge' are not allowed.", file=threshold_file)


    # 'tissue_type' data verification.
    # ******************************
    # If the optional 'tissueType' column is present, check its values are valid.
    with_tissue_type = 'tissueType' %in% col_names
    if(with_tissue_type){
        # The 'tissueType' column must be the second column of the table.
        if(col_names[2] != 'tissueType') raise_error(
            '[TissueType] must be the second column of the file.', file=threshold_file)
        # If needed, capitalize tissue type name in input data.
        input_table[,'tissueType'] = check_and_fix_tissue_type(
            tissue_type_values = input_table[,'tissueType'],
            file_values_are_from = basename(threshold_file))
    }


    # Get list of individual markers.
    # ******************************
    # marker names are retrieved from the column names of the 'thresholds_*.txt' file.
    marker_list = sort(col_names[which(! col_names %in% c('sampleName', 'tissueType'))])
    if(length(marker_list) == 0) raise_error('File has no marker columns.', file=threshold_file)

    # Check that marker names do not end in 'p' or 'm', as these are reserved symbols for
    # indicating 'positive' and 'negative' cell types.
    if(any( regexpr(pattern='p$|m$', marker_list) != -1 )) raise_error(
        "Marker names are not allowed to end in 'p' or 'm'.", file=threshold_file)


    # Create 3D array of thresholds for each sample and tissue type combination.
    # *************************************************************************
    # The dimensions of the matrix are:
    #  -> 1 = samples
    #  -> 2 = markers
    #  -> 3 = tissue types.
    threshold_table = array(data = NA,
                            dim = c(length(sample_list), length(marker_list), length(tissue_types)),
                            dimnames = list(sample_list, marker_list, tissue_types) )

    # Loop through all samples and, for each, get the threshold values.
    for(sample_name in sample_list){
        for(tissue in tissue_types){
            # Select row from the input data that matches the current sample and tissue.
            if(with_tissue_type){
                row_nb = which(input_table$sampleName == sample_name &
                                   input_table$tissueType == tissue)
            } else row_nb = which(input_table$sampleName == sample_name)

            # Verify exactly one row was selected.
            if(length(row_nb) == 0) raise_error(paste('missing threshold value for sample:',
                                                      sample_name), file=threshold_file)
            if(length(row_nb) > 1) raise_error(paste('duplicated threshold value for sample:',
                                                     sample_name), file=threshold_file)

            # Verify that all threshold values are >= 0 or == -1.
            threshold_values = as.numeric(input_table[row_nb, marker_list])
            if(any(is.na(threshold_values))){
                raise_error(paste('non-numeric threshold value for sample:', sample_name),
                            file=threshold_file)
            }
            if(any(threshold_values < 0 & threshold_values != -1)){
                raise_error(paste('threshold value < 0 and != -1 for sample:', sample_name),
                            file=threshold_file)
            }
            threshold_table[sample_name, marker_list, tissue] = threshold_values
        }


        # Check that, for a given sample, the threshold values accross tissue types are
        # either all == -1 or all != -1. In other words, we cannot have a mix of scoring
        # and phenotyping for a given marker.
        for(marker in marker_list){
            if(!(all(threshold_table[sample_name, marker,] > 0) |
                 all(threshold_table[sample_name, marker,] == -1))){
                raise_error(msg = paste0('mix of scoring and phenotyping threshold values across ',
                                         'tissues for: sample=', sample_name, ' marker=', marker),
                            file = threshold_file)
            }
        }
    }
    return(threshold_table)
}
####################################################################################################



####################################################################################################
fix_column_separator_in_table = function(input_file, write_table_to_disk = TRUE){
    # ********************************************************************************************
    # Checks that the input text file is properly formatted. The function was created to correct
    # user errors in the input "marker combination" and "thresholds" files. These input files are
    # expected to contain tab-delimited data. The function takes the following actions:
    #  - remove any leading/trailing white space on any line of the file.
    #  - check that the file does not contain a "space-tab", "tab-space" or "tab-tab" combination.
    #    If any such combination is found, it gets replaced by a single tab.
    #  - If the file contains no tabs at all, then it is assumed that one (or more) spaces have
    #    been used as delimiters instead of tabs. In this case, the first line of the file is
    #    used to estimate the number of columns in the file, and, based on this information,
    #    the spaces in the file are replaced by tabs. This makes the assumption that in the
    #    first line of the file, spaces are indeed delimiters.
    #
    # If any of the above conditions are encounterd, the content of the file is fixed and the
    # file is overwritten on disk. The function returns TRUE or FALSE to indicate success/failure.
    #
    # ********************************************************************************************
    stopifnot(file.exists(input_file))

    # Read content of input file, without any formatting expectations.
    file_connection = file(input_file, open='r', encoding=guess_file_encoding(input_file))
    file_content = readLines(con=file_connection)
    close(file_connection)
    fix_applied = FALSE
    space_tab_mix = FALSE

    # 1. Remove any line that is fully empty (or only contains spaces and/or tabs).
    offending_lines = which(file_content=='' | grepl(pattern='^[ \t]+$', x=file_content))
    if(length(offending_lines) > 0){
        cat('### WARNING: input file has empty lines [', basename(input_file), ']. \n', sep='')
        cat('###          autodesc will remove these empty lines. \n')
        # Fix content of file.
        file_content = file_content[-offending_lines]
        fix_applied = TRUE
    }
    if(length(file_content) < 2) raise_error('Input file must have >= 2 lines', file=input_file)

    # 2. Check whether file contains any leading white spaces. If yes, remove them.
    if( any(grepl(pattern='^\t|^ ', x=file_content)) ){
        cat('### WARNING: input file has one or more lines starting with tab or space [',
            basename(input_file), ']. \n', sep='')
        cat('###          autodesc will remove leading tabs and spaces. \n')
        # Fix content of file.
        file_content = sub(pattern='^[\t ]+', replacement='', x=file_content)
        fix_applied = TRUE
    }


    # 3. Check whether file contains any trailing white spaces. If yes, remove them.
    if( any(grepl(pattern='\t$| $', x=file_content)) ){
        cat('### WARNING: input file has lines ending with tab/space:',basename(input_file),'\n')
        cat('###          autodesc will remove trailing tabs and spaces. \n')
        # Fix content of file.
        file_content = sub(pattern="[\t ]+$", replacement='', x=file_content)
        fix_applied = TRUE
    }


    # 4. Check that the file does not contain a "space-tab", "tab-space" or "tab-tab" combination.
    #    If any such combination is found, it gets replaced by a single tab.
    if( any(grepl(pattern='\t\t|\t | \t|^\t|^ ', x=file_content))){
        # Fix content of file.
        file_content = gsub(pattern='[ ]+[\t]+[ ]*|[\t]+|[\t]+[ ]+[\t]*',
                            replacement='\t', x=file_content)
        fix_applied = TRUE
        space_tab_mix = TRUE
        # Note: this is an apply command that fixes the line only if it does indeed contain the
        #       tab/space mix. I ended up using the simpler option to apply the fix to all lines.
        #unlist(lapply(x, FUN=function(x) {if( grepl(pattern="\t\t|\t | \t|^\t|^ ",x=x) ){
        #                     gsub(pattern="[ ]+[\t]+[ ]*|[\t]+|[\t]+[ ]+[\t]*",
        #                     replacement='\t', x=x)
        #                     } else{x} }))

    }


    # 5. At this point, the file should contain the same number of tabs on each line, and this
    #    number should be >= 1. If this is not the case, then it probably means that the input
    #    file is badly formatted, and e.g. is using white spaces instead of tabs as delimiters.
    tab_count = unlist(lapply(gregexpr(pattern='\t', file_content),
                              FUN=function(x) sum(as.numeric(x) >= 0)))
    if( any(tab_count==0) | any(tab_count != tab_count[1]) ){

        # Replace all groups of spaces in the header with tabs. This assumes that all spaces in
        # the header are separators.
        file_content[1] = gsub(pattern='[ ]+', replacement='\t', x=file_content[1])
        tab_count[1] = sum(as.numeric(gregexpr(pattern='\t', file_content[1])[[1]]) >= 0)
        if( tab_count[1] == 0 ){
            raise_error('first line of file must contain at least 2 elements', file=input_file)
        }

        # For all lines in the text file (except for the header), we replace the last 'x'
        # sparators (tabs of spaces) with a tab. The number 'x' is the number of tab separators
        # found in the file's header.
        # We proceed as follows:
        #  -> reverse the input line.
        #  -> change all tabs to spaces (in case spaces and tabs were mixed in a same line).
        #  -> replace the first 'x' space groups with tabs: Note that the first are in fact the
        #     last, since we reversed the line.
        #  -> reverse the line again so it's back to normal.
        for(i in 2:length(file_content)){
            reversed_line = paste(rev(strsplit(file_content[i], NULL)[[1]]), collapse='')
            reversed_line = gsub(pattern='[\t]', replacement=' ', x=reversed_line)
            for(j in 1:tab_count[1]){
                reversed_line = sub(pattern='[ ]+', replacement='\t', x=reversed_line)
            }
            file_content[i] = paste(rev(strsplit(reversed_line, NULL)[[1]]), collapse='')
        }
        fix_applied = TRUE
        space_tab_mix = TRUE

        # Recompute the tab count to check that things are now OK. If the tab count is still not
        # correct, then fixing the problem is beyond autodesc's ability and an error is reported.
        tab_count = unlist(lapply(gregexpr(pattern='\t', file_content),
                                  FUN=function(x) sum(as.numeric(x) >= 0)))
        if( any(tab_count==0) | any(tab_count != tab_count[1]) ){
            raise_error(paste0('formatting error detected in input file. The file must be tab ',
                               'delimited and all rows must have the same number of elements. ',
                               'Please correct the input file and run autodesc again'),
                        file=input_file)
        }
    }

    # If something was fixed in the file content, we overwrite the existing file.
    if(fix_applied){
        # Write changes to file to disk.
        if(write_table_to_disk){
            file_connexion = file(input_file, open='wt')
            writeLines(file_content, con=file_connexion)
            close(file_connexion)
        }

        # Print warning to user.
        if(space_tab_mix){
            cat('### WARNING: unconventional column separators detected in input file [', input_file, '].\n', sep='')
            cat('###          autodesc has fixed the file, most likely successfully, but you should \n')
            cat('###          neverthless verify that the table as shown below is correct: \n')
            cat('### ******************************************************************************** \n')
            split_content = strsplit(file_content, split="\t")
            tmpDF = data.frame(matrix(unlist(split_content), ncol=length(split_content[[1]]), byrow = T),
                               row.names=NULL, stringsAsFactors = FALSE)
            colnames(tmpDF) = as.character(tmpDF[1,])
            tmpDF = tmpDF[-1,]
            print(tmpDF)
            cat('### ******************************************************************************** \n')
            cat('### \n')
        }
    }

    # Exit function and return TRUE to indicate success.
    return(TRUE)
}
####################################################################################################

