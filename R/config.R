# Global variables and function sourcing.
# Note: in principle the assignment of global variables "<<-" is not really necessary since the
# declaration is made at the top level.

# Set application version.
POSTINFORM_VERSION <<- '0.1.0'

# Define default parameter values.
DEFAULT_CELL_COMPARTMENT = 'nucleus'
DEFAULT_PHENOTYPE_CONFIDENCE_THRESHOLD = 40

# Define authorized values in input data.
NOTISSUE                 <<- 'missing'
NOTISSUES_SYNONYMS       <<- c('other', 'nothing', 'none', 'no tissue', '')
AUTHORIZED_TISSUES       <<- c('stroma', 'tumor', 'dermis', 'epidermis', 'melanocyte', 'necrosis',
                               NOTISSUE)

AUTHORIZED_COMPARTMENTS  <<- c('nucleus', 'membrane', 'cytoplasm', 'entire_cell')
AUTHORIZED_STROMA_VALUES <<- c('DAPI', 'stroma', 'other')
AUTHORIZED_TUMOR_VALUES  <<- c('CK', 'tumor')
AUTHORIZED_MARKERS       <<- c('CAL', 'CD3', 'CD4', 'CD8', 'CD11C', 'CD15', 'CD20', 'CD56', 'CD68',
                               'CD103', 'CD163', 'CD206', 'FOXP3', 'GB', 'gH2AX', 'gH2AXN', 'IDO',
                               'IL10R', 'Keratin', 'KI67', 'PD1', 'PDL1', 'PERFORIN', 'SOX10',
                               'WT1', 'CK', 'VISTA')
IGNORED_PHENOTYPES       <<- c('DAPIp', 'MISSING')

DATAREDUCE_SCRIPT      <<- file.path(dirname(dirname(sys.frame(1)$ofile)),
                                     'inst/bash/reduce_file_size.sh')
CELL_FILES_EXTENSION   <<- '_cell_seg_data.txt'
TISSUE_FILES_EXTENSION <<- '_tissue_seg_data_summary.txt'
PARAMETERS_FILE        <<- 'parameters.txt'
THRESHOLDS_FILE        <<- 'thresholds.txt'
SAMPLE_RENAME_FILE     <<- 'sample_rename.txt'
SEGMENTATION_DATA_DIR  <<- 'seg_data'
SUMMARY_OUTPUT_DIR     <<- 'summary_statistics'
LOG_FILE_NAME          <<- 'log.txt'

# Individual marker files merging thresholds.
# These parameters are the thresholds above which a warning is displayed when merging individual
# files. They have no impact on the actual output, only on the display of warning messages.
#  * MARKER_INTENSITY_THRESHOLD: difference in marker intensity values (between individual files)
#    above which a warning is shown.
#  * FILE_LOSS_PERCENTAGE_THRESHOLD: threshold (as a percentage) above which a warning is
#    displayed during file merging.
#  * SHOW_TISSUE_CATEGORY_MISMATCH_WARNING: if TRUE, the
#    "Tissue_category values differ across files" is shown. Otherwise the warning is skipped.
MARKER_INTENSITY_THRESHOLD     <<- 0.01
FILE_LOSS_PERCENTAGE_THRESHOLD <<- 5
SHOW_TISSUE_CATEGORY_MISMATCH_WARNING <<- TRUE

# InForm version data formats supported.
SUPPORTED_INFORM_VERSIONS <<- c(2.2, 2.4)

# Load source code.
source('postinform.R')
source('input_check.R')
source('load_data.R')
source('data_reduction.R')
source('individual_markers.R')
source('rename_samples.R')
source('functions.R')
source('legacy_functions.R')
