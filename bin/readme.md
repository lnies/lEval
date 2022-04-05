# BIN Modules

## Fit

## Processing

### Import files from any MR-ToF MS data format source

```
# Initialize data processor (example for MCS6Lst)
proc = MCS6Lst()
# Process list of data files
file_list = [path_to_file1, path_to_file2, ] # etc.
df = proc.process(file_list, to_csv=True) # returns dataframe with all processed .lst files. if to_csv==True, writes .csv file for each input file
# Sum all input files
df = proc.add_all(to_csv='file_name.csv') # returns 1D dataframe with the sum of all input files and a rolling sweep number. to_csv takes a file name to write summed file to a .csv file
```
Pars contains measurement parameters from .887 header file, data contains 2D array with the data, df contains dataframes saved under filename as key name

## Utilities