#!/bin/bash
####################################################################################################
# Reduce the size of Post-inForm input files.
# *******************************************
# The script loads all '*_cell_seg_data.txt' and '*_tissue_seg_data_summary.txt' files and deletes
# all unnecessary columns in them.
# This is done to increase the speed and reduce the memory footprint of the analysis in R,
# because R is slow to read large files and uses a lot of memory when opening large tables.
#
#  -> inputDir       : directory containing the data to process.
#  -> cellCompartment: cell compartement that should be considered in the analysis.
#                      e.g. 'Nucleus', 'Membrane', 'Cytoplasm', 'Entire.Cell'.
#  -> verbose        : 0 = only error message are printed. 1 = progress message are displayed.
#
# Input data check and setup.
# **************************
inputDir="$1"
outputDir="$2"
cellCompartment="$3"
verbose="$4"
[[ -z "$inputDir" ]] && echo "### ERROR: input parameter 'inputDir' is missing." && exit 1
[[ -z "$outputDir" ]] && echo "### ERROR: input parameter 'outputDir' is missing." && exit 1
[[ -z "$cellCompartment" ]] && echo "### ERROR: input parameter 'cellCompartment' is missing." \
                            && exit 1
[[ -z $verbose ]] && verbose=0
[[ $cellCompartment == "Entire.Cell" ]] && cellCompartment="Entire Cell"

# Display progress message.
if [[ $verbose -eq 1 ]]; then
	echo "################################################################################"
	echo "### Starting data reduction script."
	echo "### Input values:"
	echo "###  -> input dir       : $inputDir"
    echo "###  -> input dir       : $outputDir"
	echo "###  -> cell compartment: $cellCompartment"
	echo "### "
fi

# Create output directory if needed.
[[ "$inputDir" == "$outputDir" ]] && echo "### ERROR: 'inputDir' and 'outputDir' must differ." \
                                  && exit 1
[[ ! -d "$outputDir" ]] && mkdir -p "$outputDir"



# Reduce cell segmentation files (*_cell_seg_data.txt).
# ****************************************************
[[ $verbose -eq 1 ]] && echo "### Reduce *_cell_seg_data.txt files:"

# Retrieve all cell segmentation files.
originalIFS="$IFS"
IFS=$'\n'
fileList=( $( find "$inputDir" -type f -name '*_cell_seg_data.txt' ) )
IFS="$originalIFS"

# For each input file, keep only the columns that are needed for the analysis.
for fileName in "${fileList[@]}"; do
    [[ $verbose -eq 1 ]] && echo "###  -> ${fileName}"

    colsToKeep=$( head -n1 "$fileName" | tr '\t' '\n' | \
                  grep -i -n "^Sample Name$\|^Tissue Category$\|^Phenotype$\|^Cell ID$\|^Cell X Position$\|^Cell Y Position$\|^Confidence$\|^Annotation ID$\|^${cellCompartment} .* Mean " | \
                  cut -f1 --delimiter=':' )
    [[ -z $colsToKeep ]] && echo "### ERROR: no columns to keep found in [$fileName]" && exit 1

    # Note that in bash it's not possible to use a file as both input and output, so we create
    # a temporary file, delete the original file and then rename the temporary file.
    cut -f$( echo $colsToKeep | tr ' ' ',' ) "$fileName" > "$outputDir/$( basename $fileName )"
    unset colsToKeep
done

[[ $verbose -eq 1 ]] && echo -e "###  -> completed. \n### "
unset fileList fileName colsToKeep



# Reduce tissue surface files (*_tissue_seg_data_summary.txt).
# ***********************************************************
[[ $verbose -eq 1 ]] && echo "### Reduce *_tissue_seg_data_summary.txt files:"

# Retrieve all tissue surface files (*_tissue_seg_data_summary.txt).
IFS=$'\n'
fileList=( $( find "$inputDir" -type f -name '*_tissue_seg_data_summary.txt' ) )
IFS="$originalIFS"

# For each input file, keep only the columns that are needed for the analysis.
for fileName in "${fileList[@]}"; do
    [[ $verbose -eq 1 ]] && echo "###  -> ${fileName}"

    colsToKeep=$( head -n1 "$fileName" | tr '\t' '\n' | \
                  grep -n "^Sample Name$\|^Tissue Category$\|^Region Area\|^Annotation ID$" | \
                  cut -f1 --delimiter=':' )
    [[ -z $colsToKeep ]] && echo "### ERROR: no columns to keep found in [$fileName]" && exit 1

    cut -f$( echo $colsToKeep | tr ' ' ',' ) "$fileName" > "$outputDir/$(basename $fileName)"
    unset colsToKeep
done
[[ $verbose -eq 1 ]] && echo -e "###  -> completed. \n### "


# End of script message.
if [[ $verbose -eq 1 ]]; then
    echo "### Data reduction completed."
    echo "################################################################################"
fi
unset fileList fileName colsToKeep originalIFS
exit 0
####################################################################################################
