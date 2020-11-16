
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
#' Reclass input 'cell_values' into a binary vector indicating whether a cell is positive for a
#' given marker (1) or not (0). This function is for phenotyped markers and expects 'cell_values'
#' to be a vector of strings.
#'
#' @param marker [string] Name of marker for which to reclass cells. This must be a single marker.
#' @param cell_values [string vector] List of cell values (phenotypes) to reclassify.
#'
reclass_cells_by_marker_phenotyped <- function(marker, cell_values){

    # Verify input values.
    stopifnot(is.character(cell_values))

    # For each phenotype value, test whether maker is part of the phenotype. E.g. test whether
    # "CD8" is part of "KI67p_CD8p_CD11cp".
    marker = paste0(marker, 'p')
    return(as.numeric(sapply(strsplit(cell_values, split='_'), FUN=function(x) marker %in% x)))
}
####################################################################################################



####################################################################################################
#' Reclass input 'cell_values' into a binary vector indicating whether a cell is positive for a
#' given marker (1) or not (0). This function is for scored markers and expects 'cell_values'
#' to be a vector of numeric values, the marker intensity values for the given marker.
#'
#' @param marker [string] Name of marker for which to reclass cells. This must be a single marker.
#' @param cell_values [numeric vector] List of cell values (intensity values) to reclassify.
#' @param tissue_type [string vector] Tissue type associated to each cell. This vector must have
#'     the same length as [cell_values].
#' @param thresholds [data frame] Reclassification threshold table. Must contain a "tissue_type"
#'     column and a column bearing the name of the marker.
reclass_cells_by_marker_scored <- function(marker, cell_values, tissue_type, thresholds){

    # Verify input values.
    stopifnot(is.numeric(cell_values))
    stopifnot(length(cell_values) == length(tissue_type))
    stopifnot(sort(unique(tissue_type[tissue_type != 'missing'])) == sort(thresholds$tissue_type))

    # Reclass cell values depending on whether they are >= or < than the threshold.
    # Values >= threshold get reclassified as 1, while values < threshold get reclassified as 0.
    # The threshold value can be different for each tissue type.
    thresholds = setNames(thresholds[,marker], thresholds[,'tissue_type'])
    reclassed_values = as.numeric(sapply(1:length(cell_values),
                                         function(x) cell_values[x] >= thresholds[tissue_type[x]]))

    # If the tissue_type vector contains any "missing" values, this introduces NA values in the
    # reclassified values output. Here we set those values to 0.
    reclassed_values[which(is.na(reclassed_values))] = 0
    return(reclassed_values)
}
####################################################################################################



####################################################################################################
#' Generate a table where each row contains values for a given sample, image ID,  cell type and
#' tissue type combination. The summary_table has the following columns:
#'
#'   SampleName     ImageID TissueType             CellType Threshold SurfacePIX SurfaceMM2 CellDensity CellCount IntMean IntMedian IntMin IntMax IntSD
#'     TMA_IF02       Total     stroma                 CD8p     -1.00         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total     stroma           CD8p_total     -1.00         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total     stroma                  GBp     -1.00         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total     stroma            GBp_total     -1.00         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total     stroma                KI67p     19.99         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total     stroma          KI67p_total     19.99         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total     stroma             CD8p_GBp        NA         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total     stroma       CD8p_GBp_total        NA         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total     stroma           CD8p_KI67p        NA         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total     stroma     CD8p_KI67p_total        NA         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total     stroma                Total        NA         NA 13431416.5          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total      tumor                 CD8p     -1.00         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total      tumor           CD8p_total     -1.00         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total      tumor                  GBp     -1.00         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total      tumor            GBp_total     -1.00         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total      tumor                KI67p     19.99         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total      tumor          KI67p_total     19.99         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total      tumor             CD8p_GBp        NA         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total      tumor       CD8p_GBp_total        NA         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total      tumor           CD8p_KI67p        NA         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total      tumor     CD8p_KI67p_total        NA         NA         NA          NA        NA      NA        NA     NA     NA    NA
#'     TMA_IF02       Total      tumor                Total        NA         NA  2321735.6          NA        NA      NA        NA     NA     NA    NA
#'
#' @param sample_name [string]
#' @param image_ids [string vector]
#' @param cell_types [string vector]
#' @param tissue_list [string vector]
#' @param markers_phenotyped [string vector]
#' @param markers_scored [string vector]
#' @param thresholds [data frame]
#' @param tissue_table [data frame]
generate_summary_table <- function(sample_name, image_ids, cell_types, tissue_list,
                                   markers_phenotyped, markers_scored, thresholds, tissue_table){

    # Input check.
    stopifnot(tissue_table[,'sample_name'] == sample_name)

    # Create summary table backbone.
    # *****************************
    # Add "Total" values to the list of image IDs and cell types. These represent the total for
    # all image IDs, and the total for all cell types, i.e. all cell types together.
    image_ids  = c('Total', image_ids)
    cell_types = c(cell_types, 'Total')

    # Compute number of rows of the summary_table, i.e. image ID, cell type and tissue type
    # combinations.
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

    # Scored markers get their threshold value extracted from the thresholds table.
    for(marker in markers_scored){
    for(tissue in tissue_list){
        col_index = which(summary_table[,'CellType'] %in% paste0(marker, c('p','p_total')) &
                          summary_table[,'TissueType'] == tissue)
        summary_table[col_index, 'Threshold'] = thresholds[thresholds$tissue_type == tissue, marker]
    }
    }


    # Add surface values in mm2 to summary table.
    # ******************************************
    # Add the surface values for each image ID and tissue type to the "Total" cell type row of the
    # summary table.
    for(image_id in image_ids){
    for(tissue in tissue_list){
        row_index = which(summary_table$ImageID    == image_id &
                          summary_table$CellType   == 'Total'  &
                          summary_table$TissueType == tissue)
        surface_rows = switch(ifelse(image_id == 'Total', 1, 2),
                              which(tissue_table[,'tissue_category']==tissue),
                              which(tissue_table[,'tissue_category']==tissue &
                                    tissue_table[,'image_id']==image_id))
        summary_table[row_index, 'SurfaceMM2'] = sum(tissue_table[surface_rows,
                                                                  'region_area_surface'])
    }
    }

    # The function returns the summary_table data frame.
    return(summary_table)
}
####################################################################################################



####################################################################################################
#' Computes the following statistics for all cell types of a given image_id and tissue type:
#'
#'   - cell count      : number of cells with type "cell_type" in tissue "tissue_type".
#'   - intensity mean  : mean intensity value for cells matching "cell_type" and "tissue_type".
#'   - intensity median: median ...
#'   - intensity min   : minimum ...
#'   - intensity max   : maximum ...
#'   - intensity SD    : standard deviation ...
#'
#' The function returns a data frame with the statistic values in the same order as listed above.
#'
#' @param image_id [string] image_id for which the statistics should be computed. If this value is
#'     set to "Total" the statistics are computed on all images in cell_table.
#' @param tissue_type [string] Tissue type, e.g. "stroma", "tumor".
#' @param cell_types [string] cell type, e.g. "CD4", "CD8", "CD4_CD8".
#' @param cell_table [data frame] Data frame containing MEAN and RECLASSIFIED marker intensity
#'     values.
#'
cell_type_statistics <- function(image_id, tissue_type, cell_types, cell_table){

    # Create return data frame.
    output_df = data.frame('ImageID'    = image_id,
                           'TissueType' = tissue_type,
                           'CellType'   = c(cell_types, 'Total'),
                           'CellCount'  = 0,
                           'IntMean'    = NA,
                           'IntMedian'  = NA,
                           'IntMin'     = NA,
                           'IntMax'     = NA,
                           'IntSD'      = NA,
                           row.names    = c(cell_types, 'Total'),
                           stringsAsFactors = FALSE)

    # Subset cell_table to keep only the rows that are matching the requested tissue type and
    # image_id value. If image_id is 'Total', only subset by tissue type to get the sum across all
    # image_id values.
    # Note: it's much faster to first subset by image_id and then by tissue_type than the reverse.
    if(image_id != 'Total') cell_table = cell_table[cell_table$image_id == image_id, ]
    cell_table = cell_table[cell_table$tissue_category == tissue_type, ]

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
    for(cell_type in cell_types){

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
        row_sum_cell_type = rowSums(sub_table[, in_cell_type, drop=FALSE])
        row_sum_total     = rowSums(sub_table[, c(in_cell_type,not_in_cell_type), drop=FALSE])
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
    if(host_os == 'windows'){
        # Determine encoding by looking at the first 2 characters.
        file_connection = file(input_file, open='rb')
        file_start = paste(readBin(con=file_connection, what='raw', n=2), collapse='')
        close(file_connection)
        # UTF-8 with BOM starts with hexadecimal characters "ef bb". UTF-16 LE with "ff fe".
        if(file_start == 'efbb') return("utf-8-bom")
        if(file_start == 'fffe') return("utf-16le")
        return("utf-8")
    }
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
decompress_file <- function(input_file, allow_overwrite=FALSE, dry_run=FALSE){
    # ********************************************************************************************
    # Decompress .zip and .tar.gz input files and return the directory containing the uncompressed
    # data.
    #
    # Input parameters:
    #  -> input_file: input .tar.gz or .zip file to decompress.
    # ********************************************************************************************
    root_dir = dirname(input_file)

    # Determine archive type:
    if(endsWith(input_file, '.tar.gz')){
        type = 'tar'
    } else if(endsWith(input_file, '.zip')){
        type = 'zip'
    } else stop('Unsupported compression format. Only [.tar.gz] and [.zip] files are supported.')

    # Determine output directory:
    if(type == 'tar'){
        data_dir = file.path(root_dir, unlist(strsplit(untar(input_file, list=TRUE)[1], '/'))[1])
    } else if(type == 'zip'){
        data_dir = file.path(root_dir,
                             unlist(strsplit(zip::zip_list(zipfile=input_file)[1,1], '/'))[1])
    }
    if(dry_run) return(data_dir)

    # Test whether output already exists, and if yes delete it if allowed by user.
    if(dir.exists(data_dir)){
        if(!allow_overwrite) stop('File decompression failed: output already exists [',data_dir,']')
        unlink(data_dir, recursive=TRUE)
    }

    # Decompress file.
    if(type == 'tar'){
        untar(input_file, list=FALSE, exdir=root_dir)
    } else if(type == 'zip'){
        zip::unzip(input_file, exdir=root_dir, overwrite=TRUE)
    }

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
