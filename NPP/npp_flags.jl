# these flag variables control whether to reload the raw data from CNV files or load from saved JLD2 file.

reloadflag = 1 # set reloadflag to 1 to load from CNV, to anything other than 1 to load the saved JLD2 file.
saveflag = 1 # if the reload and saveflag are both set to 1, then the 'npp' and 'ctd' variables containing the data from the loaded CNV files are saved to a JLD2 file
