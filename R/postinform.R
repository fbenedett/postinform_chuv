####################################################################################################
#' Process a batch of samples with Post-inForm
#'
#' @description
#' Core function to run a Post-inForm analysis for a sample.
#'
#' @details
#' Core function to run a Post-inForm analysis for a sample.
#'
#' @param input_file_or_dir [string] Input file (.zip or .tar.gz) or directory to process.
#' @param command [string][optional]. One of "check", "reduce" or "all". The default is "all".
#'     -> check : only check the input data without computing statistics.
#'     -> reduce: check and reduce/standardize the input data, without computing statistics.
#'     -> all   : run the entire Post-inForm pipeline.
#' @param output_suffix [string] Suffix to use for output directory/zip file.
#' @param compress_output [logical] If TRUE, the output is compressed to a .zip file. If FALSE, the
#'     output directory is left uncompressed.
#' @param immucan_output [logical] If TRUE, the output is formatted for the IMMUCAN project. The
#'     default value is FALSE.
#' @param allow_overwrite [logical] If TRUE, existing output files or directories are silently
#'     overwritten. If FALSE, an error message is displayed if an ouput file/directory already
#'     exits.
#' @param no_bash [logical] Only applies when running the command on a Linux operating system. If
#'     TRUE, forces the data reduction process to be done in pure R (slower) rather than as a
#'     bash subprocess.
#' @return nothing.
#' @examples
#' postinform(input_file_or_dir=input_file, delete_input=FALSE, immucan_output=FALSE)
#'
postinform <- function(input_file_or_dir,
                     command = 'process',
                     output_suffix = '',
                     compress_output = TRUE,
                     immucan_output = FALSE,
                     allow_overwrite = FALSE,
                     no_bash = FALSE){

    # Argument check
    # **************
    arg_errors = c()
    # Check input file/dir.
    if(!checkmate::testString(input_file_or_dir, na.ok=F, min.chars=1, null.ok=F)) arg_errors = c(
        arg_errors, "argument 'input_file_or_dir' must be a file or directory path.")
    if(is.character(input_file_or_dir)){
        input_file_or_dir = path.expand(input_file_or_dir)
        if(!file.exists(input_file_or_dir))  arg_errors = c(
            arg_errors, paste('input file or directory could not be found:', input_file))
    }

    # Check input suffix.
    if(!checkmate::testString(output_suffix, na.ok=F, min.chars=0, null.ok=F)) arg_errors = c(
        arg_errors, "argument 'output_suffix' must be a string.")
    if(output_suffix == '') output_suffix = paste0(command,
                                                   ifelse(endsWith(command, 'e'), 'd', 'ed'))
    if(!startsWith(output_suffix, '_')) output_suffix = paste0('_', output_suffix)

    # Check input command.
    allowed_commands = c('check', 'reduce', 'process')
    if(!checkmate::testString(command, na.ok=F,
                              pattern=paste(paste0('^', allowed_commands, '$'), collapse='|'),
                              null.ok=F)) arg_errors = c(
                                arg_errors, paste0("argument 'command' must be one of: ",
                                                   paste(allowed_commands, collapse=', '), '.'))
    # Check logical input arguments.
    if(command == 'check') compress_output = FALSE
    for(x in c('compress_output', 'immucan_output', 'allow_overwrite', 'no_bash')){
        if(!checkmate::testFlag(eval(parse(text=x)), na.ok=F, null.ok=F)) arg_errors = c(
            arg_errors, paste0("argument '", x, "' must be TRUE or FALSE."))
    }

    # Check that output does not already exist. To avoid suffix duplication, remove any trailing
    # "_command" pattern from input dir. E.g. "test_reduced" becomes "test_processed" and not
    # "test_reduced_processed".
    input_dir = ifelse(file.info(input_file_or_dir)$isdir,
                       input_file_or_dir, decompress_file(input_file_or_dir, dry_run=TRUE))
    output_dir = paste0(sub('_checked$|_reduced$|_processed$', '', input_dir), output_suffix)
    output_zip = paste0(output_dir, '.zip')
    if(!allow_overwrite){
        if(file.exists(output_dir)) arg_errors = c(
            arg_errors, paste0("ouput directory already exists: ", output_dir))
        if(compress_output & file.exists(output_zip)) arg_errors = c(
            arg_errors, paste0("ouput file already exists: ", output_zip))
    }

    # Show error summary.
    if(length(arg_errors) > 0) stop("One or more error detected in input arguments:\n",
                                    paste0(' - ', arg_errors, '\n'), call.=FALSE)
    rm(arg_errors, allowed_commands)



    # Decompress input data if needed, create output directory
    # ********************************************************
    # If the input is a compressed .zip or .tar.gz file, decompress it.
    if(!file.info(input_file_or_dir)$isdir){
        message("Decompressing ", input_file_or_dir, "...")
        input_dir = decompress_file(input_file_or_dir, allow_overwrite=allow_overwrite)
    }

    # Create output directory and log file.
    if(dir.exists(output_dir)) unlink(output_dir, recursive=TRUE)
    if(!dir.create(output_dir, showWarnings=TRUE)) stop(
        paste("Unable to create output directory:", output_dir))
    LOG_FILE <<- file.path(output_dir, LOG_FILE_NAME)
    start_time = Sys.time()



    # Run the Post-inForm pipeline
    # **************************
    log_message(paste(rep('#', 80), collapse=''), padding='')
    log_message(paste0('Starting Post-inForm - version ', POSTINFORM_VERSION))
    log_message('**********************************')
    log_message(paste0('Input file     : ', input_file_or_dir))
    log_message(paste0("Output file    : ", ifelse(compress_output, output_zip, output_dir)))
    log_message(paste0('Compress output: ', ifelse(compress_output, 'yes', 'no')))
    log_message(paste0('Immucan output : ', ifelse(immucan_output, 'yes', 'no')))
    log_message(paste0('Start time     : ',
                       format(Sys.time(), '%H:%M, %d %h %Y')), add_empty_line=TRUE)

    postinform_pipeline(input_dir, output_dir, command, immucan_output, no_bash)

    log_message(paste0("Task '", command, "' completed successfully in ",
                      time_difference(start_time, Sys.time()), ' [H:M:S]'))
    log_message(paste0("Output location: ", output_dir, ifelse(compress_output, '.zip', '')))
    log_message(paste(rep('#', 80), collapse=''), padding='')


    # Final cleanup and optional data compression
    # *******************************************
    # If the original input is from a compressed file, delete the data that was extracted from it.
    if(!file.info(input_file_or_dir)$isdir) unlink(input_dir, recursive=TRUE)
    if(command == 'check') unlink(output_dir, recursive=TRUE)
    if(compress_output){
        message("Compressing and cleaning-up data...")
        if(file.exists(output_zip)) unlink(output_zip)

        # Data is compressed to .zip instead of .tar.gz to make things more windows friendly.
        wd_save = getwd()
        setwd(dirname(output_dir))
        zip::zipr(zipfile=basename(output_zip),
                  files=basename(output_dir), recurse=TRUE, compression_level=9)
        setwd(wd_save)
        rm(wd_save)

        # Delete original (uncompressed) output.
        unlink(output_dir, recursive=TRUE)
        message("Completed")
    }

    return(invisible(NULL))
}
####################################################################################################



####################################################################################################
#' Run the Post-inForm pipeline.
#'
#' @param input_dir [string] Input directory to process.
#' @param output_dir [string] Directory where to save outputs
#' @param command [string]. See description in @postinform.
#' @param immucan_output [logical] See description in @postinform.
#' @param no_bash [logical] See description in @postinform.
#' @return nothing.
#'
postinform_pipeline <- function(input_dir,
                              output_dir,
                              command,
                              immucan_output = FALSE,
                              no_bash = FALSE){

    # Input data check
    # ****************
    #  - verify that all input files needed by the application are present.
    #  - load and check marker combination file.
    #  - load and check marker thresholds file.
    #  - search for cell and tissue segmentation files in sub-directories or the session root.
    log_message('Input data check:')
    log_message(paste('input directory:', input_dir), level=2)

    #delete_unnecessary_files(input_dir)
    inputdir_check(input_dir, output_dir)
    log_message('input dir check: OK', level=2)

    input_parameters = load_session_parameters(output_dir)
    log_message('session parameters: OK', level=2)
    log_message('completed', level=2, add_empty_line=TRUE)
    if(command == 'check') return(invisible(NULL))



    # Delete unnecessary columns in segmentation files
    # ************************************************
    # The objective is to delete all columns in the input data that are not necessary for the
    # analysis (only a small subset of columns are needed). This will allow to make the the
    # loading of the files faster and consume less memory.
    log_message('Input data reduction:')
    for(subdir in list.dirs(input_dir, recursive=FALSE)){

        # Case 1: use native R code to delete unnecessary columns from the input data. This is
        # slower but works on all operating systems.
        if(guess_host_os() != 'linux' | no_bash){
            log_message('Processing files (this may take a while):')
            reduce_file_size(input_dir = subdir,
                             output_dir = file.path(output_dir, basename(subdir)),
                             cell_compartment = input_parameters$cell_compartment)

        # Case 2: the function is running on a GNU/Linux operating system. The data reduction
        # process can be done with a bash script (faster than R code).
        # Note that we need to add quotes (') around the subdir strings for the case where there
        # is a white space in the path name. This is needed here because we pass the directory
        # path string to a bash subprocess.
        } else{
            log_message(paste('processing subdirectory:', basename(subdir)), level=2)
            tmp = system2(command = DATAREDUCE_SCRIPT,
                          args = c(paste0("'", subdir, "'"),
                                   paste0("'", file.path(output_dir, basename(subdir)), "'"),
                                   input_parameters$cell_compartment,
                                   0),
                          stdout = TRUE)
            if(!is.null(attr(tmp, 'status'))) raise_error('Error in data reduction procedure.',
                                                          file=x)
        }
    }
    log_message('completed', level=2, add_empty_line=TRUE)
    if(command == 'reduce') return(invisible(NULL))



    # Input data standardization and split by sample.
    # **********************************************
    # Standardize input data so that column names are the same regardless of the version of
    # inForm. Split the data by samples.
    log_message('Input data standardisation:')
    for(subdir in list.dirs(output_dir, recursive=FALSE, full.names=TRUE)){
        for(prefix in file.path(subdir, scan_dir_for_seg_data_files(subdir))){
            log_message(paste0(basename(prefix), CELL_FILES_EXTENSION), level=2)
            standardize_and_split_cell_data(
                input_file                     = paste0(prefix, CELL_FILES_EXTENSION),
                samples                        = input_parameters$samples,
                markers_phenotyped             = input_parameters$markers_phenotyped,
                markers_scored                 = input_parameters$markers_scored,
                cell_compartment               = input_parameters$cell_compartment,
                phenotype_confidence_threshold = input_parameters$phenotype_threshold,
                delete_input_file              = TRUE)

            log_message(paste0(basename(prefix), TISSUE_FILES_EXTENSION), level=2)
            standardize_and_split_tissue_data(
                input_file        = paste0(prefix, TISSUE_FILES_EXTENSION),
                samples           = input_parameters$samples,
                delete_input_file = TRUE)
        }
        if(length(list.files(subdir)) == 0) unlink(subdir, recursive=TRUE)
    }
    log_message('completed', level=2, add_empty_line=TRUE)



    # Merge individual marker files, if any.
    # *************************************
    dirs_with_files_to_merge = sort(unique(dirname(list.files(path=file.path(output_dir),
                                                              pattern='^individualmarker_',
                                                              all.files=F, full.names=T,
                                                              recursive=T, ignore.case=F))))
    if(length(dirs_with_files_to_merge) > 0){
        log_message('Individual marker merge:')
        for(x in dirs_with_files_to_merge){
            log_message(paste0('merge data for sample ', basename(x)), level=2)
            merge_individual_marker_files(input_dir=x)
        }
        log_message('completed', level=2, add_empty_line=TRUE)
    }


    # Verify data is present for each sample.
    # **************************************
    # Verify all samples passed as input parameters are also present in the actual data.
    check_sample_names(samples_from_parameters=input_parameters$samples, output_dir)
    check_sample_directories(sample_list=input_parameters$samples, output_dir)



    # Input sample renaming.
    # *********************
    if(file.exists(file.path(output_dir, SAMPLE_RENAME_FILE))){
        log_message('Sample renaming:')

        # Modify sample names in input files.
        rename_samples(samples=input_parameters$samples, output_dir)

        # Reload input parameters to update sample names.
        input_parameters = load_session_parameters(output_dir)
        check_sample_names(samples_from_parameters=input_parameters$samples, output_dir)
        check_sample_directories(sample_list=input_parameters$samples, output_dir)

        log_message('completed', level=2, add_empty_line=TRUE)
    }



    # Progress update: print data summary report to user.
    # **************************************************
    log_message('Input data summary:')
    log_message(c(
        paste('Input directory    :', input_dir),
        paste('Number of samples  :', length(input_parameters$samples)),
        paste('                    ', input_parameters$samples),
        paste('Phenotyped markers :', paste(input_parameters$markers_phenotyped, collapse=', ')),
        paste('Scored markers     :', paste(input_parameters$markers_scored, collapse=', ')),
        paste('Tissue types       :', paste(input_parameters$tissues, collapse=', ')),
        paste('Cell compartment   :', input_parameters$cell_compartment),
        paste('Phenotype threshold:', input_parameters$phenotype_threshold),
        paste('Cell types         :', input_parameters$marker_combinations[1]),
        paste('                    ', input_parameters$marker_combinations[-1])
    ),
    level = 2, add_empty_line = TRUE)



    # Loop through all samples and process them.
    # *****************************************
    log_message('Compute cell statistics:')
    sample_output_list = list()
    for(sample_name in input_parameters$samples){
        log_message(sample_name, level=2)
        sample_output_list[[sample_name]] = run_sample(sample_name, input_parameters)
    }
    log_message('completed', level=2, add_empty_line=TRUE)



    # Save per-sample cell count and cell density summary tables to disk.
    # ******************************************************************
    log_message('Save outputs to disk:')
    log_message('Save per-sample summary files', level=2)
    for(sample_name in input_parameters$samples){
        summary_table = sample_output_list[[sample_name]]$summary_table
        summary_table = summary_table[
            which(summary_table[,'ImageID']=='Total' & summary_table[,'CellType']!='Total'),
            c('TissueType', 'CellType', 'CellCount', 'CellDensity', 'SurfaceMM2')]

        stopifnot(nrow(summary_table) == length(input_parameters$marker_combinations) *
                      length(input_parameters$tissues) * 2)
        colnames(summary_table) = c('tissue_type', 'cell_type',
                                    'cell_count', 'cell_density', 'surface_millimeters2')
        write.table(x = summary_table,
                    file = file.path(output_dir, sample_name,
                                     paste0(sample_name, '_summary_statistics.txt')),
                    sep = '\t',
                    row.names = FALSE,
                    quote = FALSE)
    }

    # Generate summary Excel file.
    if(immucan_output){
        for(sample_name in input_parameters$samples){
            # Delete tissue segmentation files. They are not needed for immucan.
            unlink(file.path(output_dir, sample_name,
                             paste0(sample_name, TISSUE_FILES_EXTENSION)))

            # Delete "sample_name" column from cell segmentation files.
            cell_seg_file = file.path(output_dir, sample_name,
                                      paste0(sample_name, CELL_FILES_EXTENSION))
            cell_seg = read.table(cell_seg_file, header=TRUE, sep='\t', as.is=T, strip.white=T)
            write.table(x=cell_seg[,-which(colnames(cell_seg) == 'sample_name')],
                        file=cell_seg_file, sep='\t', row.names=FALSE, quote=FALSE)
        }

    } else{
        log_message('generate summary spreadsheet', level=2)
        output_file_name = file.path(output_dir, paste0('summary_statistics_',
                                                             basename(output_dir),'.xlsx'))
        generate_excel_summary_file(sample_output_list, output_file_name)
    }
    log_message('completed', level=2, add_empty_line=TRUE)

    return(invisible(NULL))
}
####################################################################################################



####################################################################################################
#' Process an individual sample
#'
#' @description
#' Core function to run an analysis for a single sample. The function takes a list of
#' "input_parameters" as input, along with the name of the sample to run. The function returns
#  a "sample_object" list, whose description is given at the end of this function.
#'
#' @param sample_name [string] Name of the sample to analyze.
#' @param input_parameters [list]. List of input arguments to run the analysis.
#' @return
#'
run_sample <- function(sample_name, input_parameters){

    # Extract values from input_parameters list "object".
    # **************************************************
    session_root_dir   = input_parameters$session_root_dir
    tissue_list        = input_parameters$tissues
    markers_phenotyped = input_parameters$markers_phenotyped
    markers_scored     = input_parameters$markers_scored
    cell_types         = input_parameters$marker_combinations
    if(length(markers_scored) > 0){
        thresholds = input_parameters$thresholds
        # Keep only the thresholds for the current sample.
        thresholds = thresholds[which(thresholds[,'sample_name'] == sample_name), ]
        thresholds_subset = as.matrix(thresholds[sample_name,,])
        rownames(thresholds_subset) = dimnames(thresholds)[[2]]
        colnames(thresholds_subset) = dimnames(thresholds)[[3]]
    } else thresholds = NULL


    # Generate full list of cell types (i.e. marker combinations).
    # ***********************************************************
    # Each cell type has 2 sub-categories: the regular cell type and the "_total" cell type.
    #  - "regular" cell types:
    #    These correspond to the cells that are positive for a given marker (or combination of
    #    markers) but not for any other marker. For instance, "CD8p" corresponds to cells that
    #    are positive for CD8, but are negative for all other markers. "CD8p_CD11cp" corresponds
    #    to cells that are positive for both CD8 and CD11c, and negative for any other marker.
    #
    #  - "total" cell types:
    #    These correspond to the cells that are positive for a given marker (or combination of
    #    markers), regardless of their value for any other marker. For instance, "CD8p_total"
    #    corresponds cells that are positive for CD8, and that can be positive or negative for
    #    any other marker. In other words, "CD8p_total" corresponds to any cell that is positive
    #    for CD8. Similarly, "CD8p_CD11cp_total" corresponds to cells that are positive for both
    #    CD8 and CD11c and that can be positive or negative for any other marker.
    #
    cell_types = paste0(rep(cell_types, each=2), c('', '_total'))


    # Load cell segmentation data for current sample.
    # **********************************************
    # Load data and convert marker intensity colums to numeric values.
    cell_table = read.table(file=file.path(session_root_dir, sample_name,
                                           paste0(sample_name, CELL_FILES_EXTENSION)),
                            sep='\t', as.is=T, h=T,
                            check.names=T, strip.white=T, colClasses='character')
    col_nb = grep('cell_id|cell_[xy]_position|_mean$', colnames(cell_table))
    cell_table[, col_nb] = lapply(cell_table[, col_nb], FUN=as.numeric)
    stopifnot(!any(is.na(cell_table)))


    # Load tissue segmentation data for current sample.
    # ************************************************
    # Check that all combinations of 'sample_name', 'image_id' and 'tissue_category' in
    # surface_table are unique.
    tissue_table = read.table(file=file.path(session_root_dir, sample_name,
                                             paste0(sample_name, TISSUE_FILES_EXTENSION)),
                              sep='\t', as.is=T, h=T, strip.white=T, check.names=T,
                              colClasses=c(rep('character',3), rep('numeric',2)))
    stopifnot(ncol(tissue_table) == 5)
    check_for_duplicated_rows(tissue_table[,1:3])

    # Get the list of image ID values from the surface_table.
    image_ids = unique(tissue_table[,'image_id'])



    # Reclassify cells based on marker intensity (scoring) or phenotype data (phenotyping).
    # ************************************************************************************
    # For each marker, a cell is reclassifed as either 1 (positive for marker) or 0 (negative).
    # Reclassification is done through phenotyping for markers part of the 'markers_phenotyped'
    # list, and by scoring for markers part of the 'markers_scored' list.
    tmp = sapply(c(markers_phenotyped, markers_scored),
                 function(x) reclass_cells_by_marker(
                     marker      = x,
                     cell_values = switch(ifelse(x %in% markers_phenotyped, 1, 2),
                                          cell_table[,'phenotype'],
                                          cell_table[,paste0(x,'_mean')]),
                     thresholds  = thresholds))
    colnames(tmp) = paste(colnames(tmp), '_reclassified', sep='')
    cell_table = cbind(cell_table, tmp)



    # Generate summary table.
    # ***********************
    # Generate a table where each row corresponds to a combination of image ID, cell type and
    # tissue type.
    summary_table = generate_summary_table(sample_name, image_ids, cell_types, tissue_list,
                                           markers_phenotyped, markers_scored,
                                           thresholds, tissue_table)



    # Compute cell counts and other statistics by cell type.
    # *****************************************************
    # Get the cell count, as well as (when applicable) the marker intensity mean, median, min, max
    # and standard deviation values, for each tissue type, cell type and image ID combination.
    # Note: I timed this loop with both for() and lapply(), and for() was a few seconds faster.
    stat_table = NULL
    for(image_id in c('Total',image_ids)){
        for(tissue_type in tissue_list){
            stat_table = rbind(stat_table,
                               cell_type_statistics(cell_type=cell_types, tissue_type=tissue_type,
                                                    image_id=image_id, cell_table=cell_table))
        }
    }
    stopifnot(nrow(stat_table) == nrow(summary_table))
    stopifnot(stat_table$ImageID == summary_table$ImageID)
    stopifnot(stat_table$TissueType == summary_table$TissueType)
    stopifnot(stat_table$CellType == summary_table$CellType)
    summary_table = cbind(summary_table[,1:8], stat_table[,4:9])
    rm(stat_table)

    # Check that the cell count for '_total' cell types is >= the count of 'regular' cell types.
    for(x in which(!endsWith(summary_table[,'CellType'], '_total') &
                   summary_table[,'CellType'] != 'Total')){
        # Get the row for the matching "_total" cell type. E.g for 'CD11cp' we get the row for
        # "CD11cp_total".
        row_id = which(summary_table[,'ImageID']  == summary_table[x,'ImageID'] &
                           summary_table[,'TissueType'] == summary_table[x,'TissueType'] &
                           summary_table[,'CellType'] == paste0(summary_table[x,'CellType'], '_total'))
        # exactly one match should be found.
        stopifnot(length(row_id) == 1)
        # "_total" cell type must be >= the cell count of the regular cell type.
        stopifnot(summary_table[row_id,'CellCount'] >= summary_table[x,'CellCount'])
    }



    # Compute surface values and cell density for each cell type.
    # ***********************************************************
    # Compute surface of each cell type based on the counts of each cell type in each tissue type
    # and the total surface of each tissue type. We compute surfaces in pixels and in square
    # millimeters.
    for(image_id in c('Total', image_ids)){
        for(tissue in tissue_list){
            row_id    = which(summary_table$ImageID == image_id & summary_table$TissueType == tissue)
            row_total = which(summary_table$ImageID    == image_id &
                                  summary_table$TissueType == tissue &
                                  summary_table$CellType   == 'Total')

            # Compute surface in pixels for each cell type.
            cell_type_surface_as_proportion = switch((summary_table[row_total, 'CellCount'] == 0) + 1,
                                                     summary_table[row_id,'CellCount'] /
                                                         summary_table[row_total,'CellCount'], 0)
            total_tissue_surface = summary_table[row_total, 'SurfaceMM2']
            summary_table[row_id, 'SurfaceMM2'] = round(cell_type_surface_as_proportion *
                                                            total_tissue_surface)

            # Compute cell density for the entire tissue for each cell type.
            summary_table[row_id, 'CellDensity'] = switch((total_tissue_surface == 0) + 1,
                                                          round(summary_table[row_id,'CellCount'] /
                                                                    total_tissue_surface * 1e+06, 3), 0)
        }
    }
    # Verify no NA value is left in the surface and cell density columns of the summary table.
    stopifnot(all(!is.na(summary_table$SurfaceMM2)))
    stopifnot(all(!is.na(summary_table$CellDensity)))

    # TEMPORARY:
    # CONVERT ALL SURFACE VALUES FROM MICROMETERS TO MILIMETERS.
    summary_table[,'SurfaceMM2'] = round(summary_table[,'SurfaceMM2'] / 1e+06, 6)


    # Compute cell count statistics.
    # ******************************
    # Compute a number of statistics for the cell counts of each cell type accross the image
    # subsets of the sample. The statistics we compute are the following:
    #  - Total cell count accross all image subsets (this is also available in the summary_table).
    #  - Mean value of cell counts accross all image subsets.
    #  - Median "  "  "
    #  - Min "  "  "
    #  - Max "  "  "
    #  - Standard deviation "  "  "
    count_stat_table = summary_table[summary_table$ImageID=='Total',
                                     c('SampleName','CellType','TissueType','CellCount')]
    names(count_stat_table)[4] = 'TotalCount'

    # Create a temporary matrix that contains cell counts: rows are cell types and columns are all
    # the image subsets. we can then easily compute statistics accross the columns of that table.
    tmp_matrix = as.matrix(sapply(image_ids, FUN=function(x) summary_table[summary_table$ImageID==x,
                                                                           'CellCount']))
    count_stat_table[,'MeanCount']   = apply(tmp_matrix, MARGIN=1, FUN=mean)
    count_stat_table[,'MedianCount'] = apply(tmp_matrix, MARGIN=1, FUN=median)
    count_stat_table[,'MinCount']    = apply(tmp_matrix, MARGIN=1, FUN=min)
    count_stat_table[,'MaxCount']    = apply(tmp_matrix, MARGIN=1, FUN=max)
    count_stat_table[,'SDCount']     = apply(tmp_matrix, MARGIN=1, FUN=sd)

    # Add columns that contain the cell count total for each image subset.
    tmp_table = sapply(image_ids,
                       FUN=function(x) summary_table[summary_table$ImageID==x, 'CellCount'])
    colnames(tmp_table) = paste0('TotalCount_', colnames(tmp_table))
    stopifnot(nrow(tmp_table) == nrow(count_stat_table))

    # Verify the sum of counts accross all image subsets corresponds to the value found in the
    # "TotalCount" column.
    stopifnot(all(apply(tmp_table, MARGIN=1, FUN=sum) == count_stat_table[,'TotalCount']))

    # Add columns containing per-image total counts to the count_stat_table.
    count_stat_table = cbind(count_stat_table, tmp_table)



    # Create single-row summary tables for cell density and cell count.
    # ****************************************************************
    # The summary tables (count and density) contain a single line (the sample) and the columns
    # are the different cell types, as well as the tissue surfaces and total surface.
    summary_density = summary_count = data.frame('sampleID'=sample_name)
    col_id = 1
    tissue_cell_count   = NULL
    tissue_surface_list = NULL

    for(cell_type in c(cell_types, 'Total')){
        for(tissue in tissue_list){

            # Get row of summary_table for the current cell type.
            row_id = which(summary_table[,'ImageID']    == 'Total' &
                               summary_table[,'TissueType'] == tissue &
                               summary_table[,'CellType']   == cell_type)
            stopifnot(length(row_id) == 1)

            # Get cell density and cell count for the current cell type.
            col_id = col_id + 1
            summary_count[,col_id] = summary_table[row_id, 'CellCount']
            summary_density[,col_id] = summary_table[row_id, 'CellDensity']

            # Assign name to column.
            if(cell_type == "Total"){
                colnames(summary_count)[col_id]   = paste("total.cellCount.", tissue, sep='')
                colnames(summary_density)[col_id] = paste("total.cellDensityMM2.", tissue, sep='')

                # Update total counters for cell count and tissue surface.
                tissue_cell_count = c(tissue_cell_count, summary_table[row_id,"CellCount"])
                tissue_surface_list = c(tissue_surface_list, summary_table[row_id,"SurfaceMM2"])

            } else{
                colName = paste(sub(pattern="_total$", replacement=".total", x=cell_type), ".", tissue, sep='')
                colnames(summary_count)[col_id]   = colName
                colnames(summary_density)[col_id] = colName
                rm(colName)
            }
        }
    }

    # Add surface (in millimeters^2) for each tissue.
    col_id2 = col_id + length(tissue_list)
    col_id  = col_id + 1
    summary_count[,col_id:col_id2]   = tissue_surface_list
    summary_density[,col_id:col_id2] = tissue_surface_list
    colnames(summary_count)[col_id:col_id2]   = paste("total.SurfaceMM2.", tissue_list, sep='')
    colnames(summary_density)[col_id:col_id2] = paste("total.SurfaceMM2.", tissue_list, sep='')

    # Add total cell count and cell density accross all tissue types.
    col_id = col_id2 + 1
    summary_count[,col_id]   = sum(tissue_cell_count)
    summary_density[,col_id] = round(sum(tissue_cell_count)/sum(tissue_surface_list), 3)
    colnames(summary_density)[col_id] = "total.cellDensityMM2"
    colnames(summary_count)[col_id]   = "total.cellCount"

    # Add total tissue surface - this is the same value for both the summary_density and summary_count.
    col_id = col_id + 1
    summary_count[,col_id] = summary_density[,col_id] = round(sum(tissue_surface_list), 3)
    colnames(summary_count)[col_id] = colnames(summary_density)[col_id] = "total.surfaceMM2"



    # Create "sample" object for the current sample.
    # *********************************************
    # The sample object is a list that contains the following elements:
    #  -> sample_name     : name of the sample (a string).
    #  -> image_ids       : list of the images that constitute the sample. A sample can consist of
    #                       one or more images that were taken by the microscope for the  sample.
    #  -> summary_table   : summary table for the current sample. Contains all values for each
    #                       cell type.
    #  -> count_stat_table: statistics of cell count values across image subsets of a sample.
    #  -> summary_count   : 1-row summary for cell count. Also includes tissue surfaces.
    #  -> summary_density : single row summary for cell density. Also includes tissue surfaces.
    sample_object = list()
    sample_object$sample_name      = sample_name
    sample_object$image_ids        = image_ids
    sample_object$summary_table    = summary_table
    sample_object$count_stat_table = count_stat_table
    sample_object$summary_count    = summary_count
    sample_object$summary_density  = summary_density
    return(sample_object)
}
####################################################################################################
