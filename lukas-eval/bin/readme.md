# BIN Modules

## Fit

## Processing

### Import files from any MR-ToF MS data format source

```
# Initialize data processor (example for MCDWIN887)
proc = MCDWIN887()
# Process list of data files
file_list = [path_to_file1, path_to_file2,] # etc.
pars, data, df = proc.process(file_list, to_csv=True)
```
Pars contains measurement parameters from .887 header file, data contains 2D array with the data, df contains dataframes saved under filename as key name

## Utilities