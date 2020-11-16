# Post-inForm
Post-process cell immunofluorescence data produced by the inForm software.

**Important:** for frequently asked questions, please seet the **Frequently-made mistakes** section
on further down in this document.


##### Dependencies
The following R libraries are needed to run Post-inForm:
* openxlsx
* zip
* checkmate

These can be installed with the following R command: `install.packages(c("openxlsx", "zip", "checkmate", "stringi"))`
<br>
<br>


## Examples of how to run Post-inForm
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
* `reduce`: reduce size of input data by deleting all unnecessary data from input. For standard
  inForm data, the reduction in size is approximately of a factor 10.
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
  output file or directory are silently and mercilessly deleted. Leave this to `FALSE` to avoid
  accidental deletions.
* **output_suffix**: suffix to be appended to the input file or directory name to form the output
  name. E.g. if the input file is named `Test_session.zip` and `output_suffix` is set to
  `processed`, then the output will be named `Test_session_processed`. By default, the suffix
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
<br>
<br>


## Post-inForm input parameter file format.
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
<br>
<br>


## Frequently-made mistakes
If your post-Inform analysis fails with an error, please verify the following points:
* The input `parameters.txt` and `sample_rename.txt` files do not contain any spaces in their 
  file names names. Ideally, you will name these files simply `parameters.txt` and 
  `sample_rename.txt` but names such as `SessionName_parameters.txt` and 
  `SessionName_sample_rename.txt` are also allowed.  
  However, `SessionName Sample Rename.txt` is not a valid name because it contains spaces.
  
* The input file `parameters.txt` is based on the template file provided 
  [here](tests/parameters.txt). In particular, make sure that there are no quotation marks around 
  any line of the file.
  Here is an example of a subset of a **non-valid** `parameters.txt` file and how it should be
  corrected:
  ```
  "# List of phenotyped and scored markers."
  "phenotyped_markers: CD3, CD20, CD15, CD11c"
  
  ```
  Should be:
  ```
  # List of phenotyped and scored markers.
  phenotyped_markers: CD3, CD20, CD15, CD11c
  
  ```

* The input file `sample_rename.txt` is based on the template file provided 
  [here](tests/sample_rename.txt).
  
* The `parameter.txt` and `sample_rename.txt` files are encoded in UTF-8 format.

* If your input contains multiple files marker groups to be merged, make sure that:
    * The name of all files to be merged is the same (except for the marker names).
    * The list of markers in the file names are separated with `_`.
  Here is an example of **correct** file names:
  ```
  IMI2_Test_CD3CD15_merge_cell_seg_data
  IMI2_Test_CKCD11c_merge_cell_seg_data
  ```
  And here are examples of **wrong** file names:
    * In the first file name, the "CD3CD15" marker names are not properly separated with 
      `_` characters.
    * In the second file name, a "_2" was added to the second file name that is not present in the 
      first file name. Therefore post-inForm will not be able to identify that these two files 
      should be merged.
  ```
  IMI2_TestCD3CD15_merge_cell_seg_data
  IMI2_Test_2_CKCD11c_merge_cell_seg_data
  ```
  
* `Split by coordinate failed` error: getting this error means that the splitting + renaming of the 
  specified sample failed because it could not be automatically split (vertically) in two distinct 
  samples. Please check the following:
    * The sample really contains 2 samples that must be split. Very often, this error occurs 
      because a sample that is indicated as having to be split (i.e. it has 2 `new_name` values 
      in the `sample_rename.txt` file), is in fact only a single sample and should not be split. 
      **Solution:** edit the `sample_rename.txt` file so that only one "new name" value is present.

  
