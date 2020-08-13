
####################################################################################################
check_for_duplicated_rows <- function(data_frame){
    # ********************************************************************************************
    # Checks for duplicated rows in the input data frame and raises and error if duplicates are
    # found.
    # ********************************************************************************************
    duplicated_rows = which(duplicated(data_frame, MARGIN=1))
    if(length(duplicated_rows) > 0) raise_error(
        msg = "Duplicated rows found in input data:",
        items_to_list = sapply(duplicated_rows,
                               FUN=function(x) paste(data_frame[x,], collapse=' ')))
}
####################################################################################################



####################################################################################################
reclass_cells_by_marker <- function(marker,
                                    cell_values,
                                    thresholds = NULL){
    # ********************************************************************************************
    # Computes a binary vector indicating whether a cell is positive for a given marker (1) or
    # not (0). The function accepts cell_values that are either phenotypes, or marker intensity
    # values.
    #
    # Input arguments:
    #   marker     : string. Name of marker for which to reclassify cells.
    #   cell_values: if the marker was phenotyped: must be a vector of strings.
    #                if the marker was scored: must be a vector of floats.
    #   thresholds : list of thresholds for the given marker. If the marker was phenotyped (as
    #                opposed to scored), then thresholds must be set to NULL. This is how the
    #                function knows that the maker is phenotyped.
    # ********************************************************************************************

    # Case 1: marker was phenotyped.
    # *****************************
    if(length(thresholds) == 0){
    # Test whether the elements of the cell_values vector contain the string "markerName + p".
    # E.g. for CD11c, we test if the phenotype contains the string "CD11cp".
    # "partial matches", means that if we are searching for say "CD11", then all cells
    #                    whose phenotype contains "CD11p" (e.g. "CD11p", "CD11p_CD86p" or "CD8p_CD11p_CD86p") will be
    #                    reclassified as 1 (the marker is present in the cell). Note that if the marker apprears in the
    #                    phenotype with a "m", e.g. "CD8p_CD11m", then this is not counted as a match since the "m"
    #                    stands for "minus", meaning the marker is not present in the cell.
        marker = paste0(marker, 'p')
        stopifnot(is.character(cell_values))
        return(sapply(strsplit(cell_values, split='_'), FUN=function(x) ifelse(marker %in% x,1,0)))
    }


    # Case 2: marker was scored.
    # **************************
    # In this case we reclass the values depending on whether they are >= or < than the threshold.
    # Values >= threshold get reclassified as 1, while values < threshold get reclassified as 0.
    stopifnot(is.numeric(cell_values))

    stop("########### NOT IMPLEMENTED YET !!!!")
    # Create temporary vector of thresholds. This vector has the same length as the number of rows in
    # the cell table, and simply contains the threshold value coresponding to the tissue type for
    # each row of the cell table.
    thresholdVector = as.numeric(sapply(intensity_table[,'tissue_category'],
                                        FUN=function(x) thresholdList[1,x]))

    # Test whether intensity value is >= threshold.
    return(as.numeric(intensity_table[,tmpColNb] > thresholdVector))
}
####################################################################################################



####################################################################################################
generate_summary_table <- function(sample_name,
                                   image_ids,
                                   cell_types,
                                   tissue_list,
                                   markers_phenotyped,
                                   markers_scored,
                                   thresholds,
                                   tissue_table){
    # ********************************************************************************************
    # Generate a table where each row contains values for a given sample, image ID,  cell type and
    # tissue type combination. The summary_table has the following columns:
    #  - sample_name : name of sample.
    #  - image_id    : name of image subset within the sample.
    #  - tissue_type : name of tissue (e.g. "stroma", "tumor").
    #  - cell_type    : name of cell type. Either an individual marker or a combination (e.g. "CD4", "CD8p_FPXP3n").
    #  - threshold   :  threshold value for mean marker intensity values. These are the threshold that are used to
    #                  reclassify the marker intensity values into 0/1 values. Note that these values only exist
    #                  for the rows corresponding to individual markers (combination cell types and Total rows
    #                  have no values in the column).
    #  - SurfacePIX : surface of each cellType in pixel units.
    #  - SurfaceMM2 : surface of each cellType, in square mm.
    #  - CellDensity: density of cells per square mm.
    #  - CellCount  : cell count for the given cell type in the given tissue type.
    #  - IntMean    : mean value of marker intensity for the cell belonging to the row's cellType.
    #  - IntMedian  : median value   "  "  "
    #  - IntMin     : minimum value  "  "  "
    #  - IntMax     : maximum value  "  "  "
    #  - IntSD      : standard deviation value "  "  "
    #
    #   SampleName      ImageID     CellType TissueType Threshold SurfacePIX  SurfaceMM2
    #    His2757_7        Total          CD4     Stroma       0.5         NA          NA
    #    His2757_7        Total          CD4      Tumor       0.5         NA          NA
    #    His2757_7        Total        FOXP3     Stroma       0.2         NA          NA
    #    His2757_7        Total        FOXP3      Tumor       0.2         NA          NA
    #    His2757_7        Total  CD4p_FOXP3p     Stroma        NA         NA          NA
    #    His2757_7        Total  CD4p_FOXP3p      Tumor        NA         NA          NA
    #    His2757_7        Total        Total     Stroma        NA     299488          NA
    #    His2757_7        Total        Total      Tumor        NA    2108192          NA
    #    His2757_7  39995_14773          CD4     Stroma       0.5         NA          NA
    #    His2757_7  39995_14773          CD4      Tumor       0.5         NA          NA
    #    His2757_7  39995_14773        FOXP3     Stroma       0.2         NA          NA
    #    His2757_7  39995_14773        FOXP3      Tumor       0.2         NA          NA
    #    His2757_7  39995_14773  CD4p_FOXP3p     Stroma        NA         NA          NA
    #    His2757_7  39995_14773  CD4p_FOXP3p      Tumor        NA         NA          NA
    #    His2757_7  39995_14773        Total     Stroma        NA     149744          NA
    #    His2757_7  39995_14773        Total      Tumor        NA    1054096          NA
    #
    #
    # Input arguments:
    #   sample_name: string. Name of sample.
    #
    # ********************************************************************************************
    stopifnot(all(tissue_table[,'sample_name'] == sample_name))

    # Create summary table backbone.
    # *****************************
    # Add "Total" values to the list of image IDs and cell types. These represent the total for
    # all image IDs, and the total for all cell types (i.e. all cell types together).
    image_ids  = c('Total', image_ids)
    cell_types = c(cell_types, 'Total')

    # Compute nuber of rows of the summary_table (i.e. image ID, cell type and tissue type combinations)
    row_nb = length(tissue_list) * length(cell_types) * length(image_ids)

    # Create table.
    summary_table = data.frame(
            SampleName  = rep(sample_name, times=row_nb),
            ImageID     = rep(image_ids, each=length(cell_types) * length(tissue_list)),
            TissueType  = rep(sort(tissue_list), each=length(cell_types), times=length(image_ids)),
            CellType    = rep(cell_types, each=1, times=length(tissue_list) * length(image_ids)),
            Threshold   = NA,
            SurfacePIX  = NA,
            SurfaceMM2  = NA,
            CellDensity = NA,
            CellCount   = NA,
            IntMean     = NA,
            IntMedian   = NA,
            IntMin      = NA,
            IntMax      = NA,
            IntSD       = NA,
            stringsAsFactors = FALSE)


    # Add threshold values to summary_table.
    # *************************************
    # Threshold values for phenotyped markers are set to -1
    summary_table[summary_table[,'CellType'] %in%
                  paste0(rep(markers_phenotyped, each=2), c('p','p_total')), 'Threshold'] = -1

    # Scored markers get the threshold value extracted from the thresholds table.
    for(marker in markers_scored){
    for(tissue in tissue_list){
        stop("######## NOT IMPLEMENTED YET")
        tmpColNb = which( summary_table[,'CellType'] %in% paste0(marker, c('p','p_total'), sep='') &
                          summary_table[,'TissueType'] == tissue )
        summary_table[tmpColNb,'Threshold'] = threshold_table[sample_name,marker,tissue]
    }
    }


    # Add surface values in MM2 to summary table.
    # ******************************************
    # Add the surface values for each image ID and tissue type to the "Total" cell type row of the
    # summary table.
    for(image_id in image_ids){
    for(tissue in tissue_list){
        row_nb = which(summary_table$ImageID    == image_id &
                       summary_table$CellType   == 'Total'  &
                       summary_table$TissueType == tissue)
        surface_rows = switch(ifelse(image_id == 'Total', 1, 2),
                              which(tissue_table[,'tissue_category']==tissue),
                              which(tissue_table[,'tissue_category']==tissue &
                                    tissue_table[,'image_id']==image_id))
        summary_table[row_nb, 'SurfaceMM2'] = sum(tissue_table[surface_rows, 'region_area_surface'])
    }
    }

    # The function returns the summary_table data frame.
    return(summary_table)
}
####################################################################################################



####################################################################################################
cell_type_statistics <- function(image_id, tissue_type, cell_types, cell_table){
    # *******************************************************************************************
    # Computes the following statistics for all cell types of a given image_id and tissue type:
    #  - cell count      : number of cells with type "cell_type" in tissue "tissue_type".
    #  - intensity mean  : mean intensity value for cells maching "cell_type" and "tissue_type".
    #  - intensity median: median ...
    #  - intensity min   : minimum ...
    #  - intensity max   : maximum ...
    #  - intensity SD    : standard deviation ...
    #
    # The function returns a data frame with the statistic values in the same order as listed above.
    #
    # Input arguments:
    #   cell_type  : string describing cell type, e.g. "CD4", "CD8", "CD4_CD8".
    #   tissue_type: string for tissue type, e.g. "stroma", "tumor".
    #   cell_table : data frame containing MEAN and RECLASSIFIED marker intensity values,
    #                      in the format as returned by the load_cell_data() function.
    #   image_id   : string giving the image_id value of the subset image for which the
    #                statistics should be computed. If this value is set to "Total" (the
    #                default), then the statistics are computed on all images contained in
    #                the cell_table.
    # *******************************************************************************************

    # Create return data frame.
    cell_types = c(cell_types, 'Total')
    output_df = data.frame('ImageID'    = image_id,
                           'TissueType' = tissue_type,
                           'CellType'   = cell_types,
                           'CellCount'  = 0,
                           row.names    = cell_types,
                           stringsAsFactors = FALSE)
    output_df[,c('IntMean', 'IntMedian', 'IntMin', 'IntMax', 'IntSD')] = NA

    # Subset cell_table to keep only the rows that are matching the requested "tissue_category"
    # and "image_id" values. If the image ID is 'Total', then we only subset by tissue type
    # because it corresponds values for the sum accross all image_id values.
    # Note: it's much faster to first subset by image_id and then by tissue_type than the reverse.
    # get("image_id", 1L)
    if(image_id != 'Total') cell_table = cell_table[cell_table$image_id==image_id, ]
    cell_table = cell_table[cell_table$tissue_category==tissue_type, ]
    # If the subset has 0 rows, the image ID has no cells for the given tissue type. This does
    # sometimes occur.
    if(nrow(cell_table) == 0) return(output_df)
    output_df['Total', 'CellCount'] = nrow(cell_table)

    # Tests with data.table.
    #img_id = image_id
    #if(img_id != 'Total') sub_table = cell_table[image_id == img_id]
    #sub_table = sub_table[tissue_category == tissue_type]
    #
    #system.time(dt[image_id == tmp])
    #system.time(cell_table[cell_table$image_id==image_id, ])
    #system.time(dt[tissue_category == tissue_type])
    #system.time(cell_table[cell_table$tissue_category==tissue_type, ])


    # Loop through all cell types and compute statistics
    for(cell_type in cell_types[-length(cell_types)]){

        # Split the cell type into individual markers.
        # *******************************************
        # cell_type are strings that contain the names of individual markers, followed by either
        # 'p' (positive) or 'm' (negative), separated by '_' characters, e.g: "CD68p_CD11cp"
        tmp = unlist(strsplit(sub('_total$', '', cell_type), split='_'))
        marker_names    = as.character(sapply(tmp, FUN=function(x) substr(x, 1, nchar(x)-1)))
        marker_statuses = as.character(sapply(tmp, FUN=function(x) substr(x, start=nchar(x),
                                                                          stop=nchar(x))))
        stopifnot(all(marker_statuses %in% c('p','m')))

        # For any marker with negative status ("m"), inverse the values in the reclassified
        # column of the sub_table: the 0 become 1 and 1 become 0.
        sub_table = cell_table[]
        for(x in marker_names[marker_statuses == 'm']) sub_table[, paste0(x,'_reclassified')] =
            1 - sub_table[, paste0(x,'_reclassified')]

        # Identify which reclassified columns are part of the cell type or not.
        col_names = colnames(sub_table)
        in_cell_type = which(col_names %in% paste0(marker_names, '_reclassified'))
        not_in_cell_type = which(endsWith(col_names, '_reclassified') &
                                 !col_names %in% col_names[in_cell_type] )


        # Subset intensity values table to keep only the rows that are matching the cell_type.
        # ***********************************************************************************
        # The way to do this depends on whether the cell type is a "regular" or '_total' type:
        #  -> total  : keep all cells (rows of table) that are positive for the markers associated
        #              to the target cell_type.
        #  -> regular: keep only cells (rows of table) that are positive for the markers of the
        #              target cell_type and negative for all other markers.
        row_sum_cell_type = rowSums(sub_table[, in_cell_type, drop=F])
        row_sum_total     = rowSums(sub_table[, c(in_cell_type,not_in_cell_type), drop=F])
        if(endsWith(cell_type, '_total')){
            sub_table = sub_table[row_sum_cell_type == length(in_cell_type),]
        } else{
            sub_table = sub_table[row_sum_cell_type == length(in_cell_type) &
                                  row_sum_cell_type == row_sum_total, ]
        }


        # Compute cell count and intensity value statistics.
        # *************************************************
        # Cell count.
        output_df[cell_type, 'CellCount'] = nrow(sub_table)

        # Intensity value statistics are only computed for single marker types. In addition,
        # phenotyped markers don't necessarily have an associated intensity column.
        col_nb = which(col_names %in% paste0(marker_names, '_mean'))
        if(length(marker_names) == 1 & length(col_nb) == 1){
            if(nrow(sub_table) > 0){
                output_df[cell_type, 'IntMean']   = round(mean(sub_table[,col_nb]), 3)
                output_df[cell_type, 'IntMedian'] = round(median(sub_table[,col_nb]), 3)
                output_df[cell_type, 'IntMin']    = round(min(sub_table[,col_nb]), 3)
                output_df[cell_type, 'IntMax']    = round(max(sub_table[,col_nb]), 3)
                output_df[cell_type, 'IntSD']     = round(sd(sub_table[,col_nb]), 3)
            } else output_df[cell_type, 5:9] = 0
        }
    }
    return(output_df)
}
####################################################################################################



####################################################################################################
guess_host_os = function(){
    # ********************************************************************************************
    # Determine what is the operating system of the machine executing this function.
    #
    # ********************************************************************************************
    # Get name of host operating system.
    sysinf = Sys.info()

    # Case 1: name of opperating system can be retrieved with "Sys.info()". This is the usual case.
    if(!is.null(sysinf)){
        os <- sysinf['sysname']
        if(os == 'Darwin'){
            os = 'osx'
        } else if(os == 'Linux'){
            os = 'linux'
        } else if(os == 'Windows'){
            os = 'windows'
        } else os = 'unknown'

    # Case 2: Exotic OS and Sys.info() did not work.
    } else{
        os = .Platform$OS.type
        if(grepl('^darwin', R.version$os)){
            os = 'osx'
        } else if(grepl('linux-gnu', R.version$os)){
            os = 'linux'
        } else os = 'unknown'
    }
    return(os)
}
####################################################################################################



####################################################################################################
#' Determine encoding of input file
#'
#' @param input_file
#' @return [string] file encoding.
guess_file_encoding = function(input_file){

    # Make a system call to identify the encoding of input_file. Since the system call depends on
    # the host OS, so we have to determine that first.
    host_os = guess_host_os()
    if(host_os == 'linux') tmp = system2("file", args=c("-i", paste0("'",input_file,"'")), stdout=T)
    if(host_os == 'osx') tmp = system2("file", args=c("-I", paste0("'",input_file,"'")), stdout=T)
    if(host_os == 'windows') return("native.enc")
    if(host_os == 'unknown') raise_error("Unable to detect host OS.")
    encoding_type = sub(pattern='^.* charset=', replacement='', x=tmp)

    # Verify that the encoding is part of a pre-determined list.
    encoding_list = c('us-ascii', 'ascii', 'utf-8', 'utf-16le')
    if(! encoding_type %in% encoding_list) raise_error(msg='Unable to identify encoding of file.',
                                                       file=input_file)
    return(encoding_type)
}
####################################################################################################



####################################################################################################
generate_global_summary_table = function(sample_output_list,
                                         cell_measurement,
                                         column_order){
    # ********************************************************************************************
    # Generate a global summary table for cell counts or cell densities for all samples.
    # The user can select to have columns ordered by cell type or by tissue type.
    #
    # Input arguments:
    #   sample_output_list: list of "sample output" objects, as produced by "run_sample" function.
    #   cell_measurement  : either 'count' or 'density'.
    #   column_order      : either 'by_cell_type' or 'by_tissue'.
    # ********************************************************************************************
    stopifnot(length(sample_output_list) > 0)
    required_items = c('summary_count','summary_density', 'summary_table')
    stopifnot(sapply(sample_output_list, FUN=function(x) all(required_items %in% names(x)) ))


    # Generate summary tables for cell counts/density.
    # ***********************************************
    # The table has one row for each sample. The columns contain cell counts or cell density
    # values for all cell type and tissue type combinations.
    table_name = ifelse(cell_measurement == 'count', 'summary_count', 'summary_density')
    summary_table = do.call(rbind.data.frame,
                            lapply(sample_output_list, FUN=function(x) x[[table_name]]))

    if(column_order == 'by_cell_type') return(summary_table)


    # If requested, order columns by tissue type.
    # ******************************************
    total_col = ncol(summary_table)
    tissue_count = length(unique(sample_output_list[[1]]$summary_table[,'TissueType']))
    col_order = c(1,
                  unlist(lapply(1:tissue_count,
                                FUN=function(x) seq(from=1+x, to=total_col-2, by=tissue_count))),
                  total_col - 1,
                  total_col)
    stopifnot(length(col_order) == total_col)
    return(summary_table[,col_order])
}
####################################################################################################



####################################################################################################
generate_excel_summary_file = function(sample_output_list, output_file_name){
    # ********************************************************************************************
    # Create a summary file in Microsoft Excel format.
    # The excel file contains the 'summary' and 'count' tables for each sample in separate
    # worksheets. It also contains global summary tables for cell counts and densities.
    #
    # ********************************************************************************************
    stopifnot(length(sample_output_list) > 0)

    # Add global summary tables to the Excel file (i.e. workbook).
    # ***********************************************************
    # Compute global (i.e. for all samples) cell count and cell density tables. Add them to the
    # Excel workbook.
    wb = openxlsx::createWorkbook()
    for(cell_measurement in c('count', 'density')){
        for(column_order in c('by_cell_type', 'by_tissue')){

            worksheet_name = paste0('summary_', cell_measurement)
            if(column_order == 'by_tissue') worksheet_name = paste0(worksheet_name, '_v2')

            openxlsx::addWorksheet(wb=wb, sheetName=worksheet_name)
            openxlsx::writeData(wb=wb, sheet=worksheet_name, borders='none',
                                x=generate_global_summary_table(sample_output_list,
                                                                cell_measurement, column_order))
        }
    }

    # Add per sample summary and cell count tables to workbook.
    # ********************************************************
    # Define name of worksheet. A complication is that Excel worksheet names cannot exeed 31
    # characters, so depending on the maximum length of sample names, the worsheet name will be
    # named either:
    #  - sample name + '_summary'/'_count'
    #  - 'sample' + incremental number + '_summary'/'_count'   (e.g. sample1_summary).
    prefix_list = names(sample_output_list)
    if(any(nchar(prefix_list) > 23)) prefix_list = paste0('sample', 1:length(prefix_list))

    for(x in 1:length(sample_output_list)){
        for(table_name in c('summary_table','count_stat_table')){

            worksheet_name = paste0(prefix_list[x], '_', sub('_.*$', '', table_name))
            openxlsx::addWorksheet(wb=wb, sheetName=worksheet_name)
            openxlsx::writeData(wb=wb, sheet=worksheet_name,
                                x=sample_output_list[[x]][[table_name]], borders='none')
        }
    }

    # Write Excel workbook to disk.
    openxlsx::saveWorkbook(wb, output_file_name, overwrite=TRUE)
}
####################################################################################################



####################################################################################################
verify_columns_are_present = function(required_columns, actual_columns, file_name){
    # ********************************************************************************************
    # Verifies that all values listed in required_columns are present in actual_columns. If not
    # then an error is raised.
    #
    # Input arguments:
    #  - required_columns: string vector. List of column names that are expected to be present.
    #  - actual_columns: string vector. List of column names to verify.
    #  - file_name: name of file from where the column names are originally taken. This is only
    #               used to display the file name in the error message.
    # ********************************************************************************************
    missing_columns = which(! required_columns %in% actual_columns)
    if(length(missing_columns) > 0) raise_error(
        msg = 'The following columns are missing in the input file:',
        file = file_name,
        items_to_list = required_columns[missing_columns])
}
####################################################################################################



####################################################################################################
extract_imageid <- function(input_vector, input_file='file not specified'){
    # ********************************************************************************************
    # Extract "image ID" values from input_vector strings that have the following structure
    # depending on the version of inForm:
    # inForm 2.2. "image ID" values are found in the "Sample Name" column.
    #  - <sample name>_[<image ID>].im3 , e.g. "His2757_7_[39995,14773].im3".
    # inForm 2.4. "image ID" values are found in the "Annotation ID" column.
    #  - <sample name>_<some text>_[<image ID>], e.g. "Bern_TMA1_Scan1_Core[1,F,1]_[11914,31856]".
    #
    # Note that sample names are not allowed to to contain '[' or ']' characters, because these
    # are reserved to be used as delimiters of the image_id.
    #
    # Input parameters:
    #  -> input_vector: vector of strings from which to extract sample_name or image_id values.
    #  -> input_file:   name of file from which the input strings were loaded from. This is only
    #                   needed to display an error message.
    # ********************************************************************************************
    stopifnot(length(input_vector) > 0)

    # Remove trailing ".im3" or "_Mx.im3" suffix found in sample name values of inForm 2.2 files.
    input_vector = sub('(_M[0-9])?.im3$', '', input_vector)

    # Verify that the input strings contain at least one '[' and one ']' character. This is
    # mandatory because the general structure of the strings must be '<sample name>_[xxxxx,xxxxx]'.
    if(!all(grepl(pattern='_\\[', input_vector) & endsWith(x=input_vector, suffix=']'))){
        raise_error(msg = "One or more values do not follow the format <sample name>_[xxxxx,xxxxx]",
                    file = input_file)
    }

    # Keep only the image ID part of the string and replace the ',' separator with '_'.
    input_vector = sub(pattern='^.*_\\[', replacement='', x=input_vector)
    input_vector = sub(pattern=']$', replacement='', x=input_vector)
    input_vector = sub(pattern=',', replacement='_', x=input_vector)
    return(input_vector)
}
####################################################################################################



####################################################################################################
decompress_file <- function(input_file, dry_run=FALSE){
    # ********************************************************************************************
    # Decompress .zip and .tar.gz input files and return the directory containing the uncompressed
    # data.
    #
    # Input parameters:
    #  -> input_file: input .tar.gz or .zip file to uncompress.
    # ********************************************************************************************
    root_dir = dirname(input_file)

    # Input file is tarball.
    if(endsWith(input_file, '.tar.gz')){
        data_dir = file.path(root_dir, unlist(strsplit(untar(input_file, list=TRUE)[1], '/'))[1])
        if(dry_run) return(data_dir)
        if(dir.exists(data_dir)) unlink(data_dir, recursive=TRUE)
        untar(input_file, list=FALSE, exdir=root_dir)

    # Input file is zip archive.
    } else if(endsWith(input_file, '.zip')){
        data_dir = file.path(root_dir,
                             unlist(strsplit(zip::zip_list(zipfile=input_file)[1,1], '/'))[1])
        if(dry_run) return(data_dir)
        if(dir.exists(data_dir)) unlink(data_dir, recursive=TRUE)
        zip::unzip(input_file, exdir=root_dir, overwrite=TRUE)

    # Unsupported compression format.
    } else stop('Unsupported compression format.')


    # Return dirname of extracted files.
    stopifnot(file.info(data_dir)$isdir)
    return(data_dir)
}
####################################################################################################



####################################################################################################
time_difference <- function(start, end){
    # ********************************************************************************************
    #
    # ********************************************************************************************
    diff_in_seconds = round(as.numeric(end) - as.numeric(start))
    hours = floor(diff_in_seconds / 3600)
    minutes = floor((diff_in_seconds - (hours * 3600)) / 60)
    seconds = diff_in_seconds - (hours * 3600) - (minutes * 60)
    minutes = ifelse(nchar(minutes) == 1, paste0('0', minutes), minutes)
    seconds = ifelse(nchar(seconds) == 1, paste0('0', seconds), seconds)
    return(paste0(hours, ':', minutes, ':', seconds))
}
####################################################################################################



####################################################################################################
log_message <- function(message, level=1, padding='### ',
                        add_empty_line = FALSE, log_to_file=TRUE, print_to_stdout=TRUE){
    # ********************************************************************************************
    #
    # ********************************************************************************************
    message = paste0(padding, switch(level,'','-> '), message)
    if(add_empty_line) message = c(message, padding)

    if(log_to_file) write(message, LOG_FILE, append=T)
    if(print_to_stdout) cat(message, sep='\n')
}
####################################################################################################



####################################################################################################
raise_error <- function(msg, file='', items_to_list=NULL, type='error'){
    # ********************************************************************************************
    # Displays am error or warning message.
    #
    # Input arguments:
    #  msg          : string vector. One or more error messages to display.
    #  file         : string. Optional file in which the error occurred.
    #  items_to_list: string vector. Items to list in the error message.
    #  type         : either 'error' or 'warning'.
    # ********************************************************************************************
    if(type == 'error'){
        prefix = 'ERROR: '
        log_padding = ''
    } else{
        prefix = 'WARNING: '
        log_padding = '### '
    }
    padding = strrep(' ', nchar(prefix))

    # Display name of file where problem was detected, if applicable.
    if(file != '') log_message(paste0(prefix, 'problem detected in ',
                                      ifelse(dir.exists(file), 'directory', 'file'),
                                      ': ', basename(file)), padding = log_padding)

    # Display error/warning message.
    log_message(paste0(c(ifelse(file=='', prefix, padding), rep(padding, length(msg)-1)), msg),
                padding = log_padding)

    # List items, if any.
    if(length(items_to_list) > 0) log_message(paste0(padding, ' - ', items_to_list),
                                              padding = log_padding)

    if(type == 'error')  stop('error detected. See message above.', call.=FALSE)
}
####################################################################################################
