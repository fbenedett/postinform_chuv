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
