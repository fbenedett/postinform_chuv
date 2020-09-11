# Post-inForm
Post-process cell immunofluorescence data produced by the inForm software.

##### Dependencies
The following R libraries are needed to run Post-inForm:
* openxlsx
* zip
* checkmate

These can be installed with the following R command: `install.packages(c("openxlsx", "zip", "checkmate"))`


### Examples of how to run Post-inForm
##### Load Post-inForm source code:
Download or clone the project's directory to your local machine. Then set `POSTINFORM_ROOT` to the
`postinform` directory on your local machine and run the `source()` command in your R console, as
shown here:
```
POSTINFORM_ROOT = '~/projects/chuv_cte/postinform.git'
source(file.path(POSTINFORM_ROOT, 'R', 'config.R'), chdir=TRUE)
```


##### Process samples:
Post-inForm can be run with 3 different commands:
* `check`: check input data only. Does not produce any output.
* `reduce`: reduce size of input data by deleting all unecessary data from input. For standard
  inForm data, the reduction in size is approximatively of a factor 10.
* `process`: run the full Post-inForm data processing pipeline.

```
input_file = 'Session_Batch6_1_Curry2.zip'
postinform(input_file_or_dir=input_file, command='check')
postinform(input_file_or_dir=input_file, command='reduce')
postinform(input_file_or_dir=input_file, command='process')
```


##### Post-inForm options:
* **compress_output**: if `TRUE`, the output is compressed to a .zip file. If `FALSE` the output
  remains in an uncompressed directory.
* **allow_overwrite**: if `TRUE`, pre-existing files and directories with the same name as an
  output file or directory are silently and mercylessly deleted. Leave this to `FALSE` to avoid
  accidental deletions.
* **output_suffix**: suffix to be appended to the input file or directory name to form the output
  name. E.g. if the input file is named `Test_session.zip` and `output_suffix` is set to
  `processed`, then the ouput will be named `Test_session_processed`. By default, the suffix
  value is set to `reduced` when running the "reduce" command, and `processed` when running the
  "process" command.
* **immucan_output**: if `TRUE`, produces IMMUCAN compatible outputs.

Examples:  
This command will produce an uncompressed output directory "Test_session_reduced". If such
a directory already exists, it will be overwritten since "allow_overwrite" is set to TRUE
```
postinform(input_file_or_dir="Test_session.zip", command='reduce',
         compress_output=F, immucan_output=FALSE, allow_overwrite=TRUE)
```

This command will produce an output file named "Test_session_random_suffix.zip".
```
postinform(input_file_or_dir="Test_session.zip", command='process', output_suffix="random_suffix",
         compress_output=TRUE, immucan_output=TRUE, allow_overwrite=FALSE)
```



### Post-inForm input parameter file format.
The list of samples, tissues, markers and marker combinations to process are passed to post-inForm 
via a single plain text file that must be named `parameters.txt` and be located at the root of the 
input directory.
**Important:** for **Windows users**, the `parameters.txt` file should be encoded in `UTF-8` 

The parameters input file must contain the following 5 sections. For each section, values can be 
passed on the same line separated by a `,`, or on multiple lines (one value per line). Any line 
starting the a `#` value is ignored (allows to add comments to file).  
For the `marker_combinations:` section, the special value `all` can be passed to process all 
possible combinations of markers (to avoid having to enter all combinations manually).
```
samples:
tissues:
phenotyped_markers:
scored_markers: 
marker_combinations: 
```

A template file can be downloaded [here](tests/parameters.txt).  
`parameters.txt` file example:
```
# List of samples to process.
samples:
SAMPLE_1
SAMPLE_2
SAMPLE_3

# List of tissues to process.
tissues: stroma, tumor

# List of markers.
phenotyped_markers: CD20, CD3, CD68
scored_markers: 

# Marker combinations to test
marker_combinations: all
```



