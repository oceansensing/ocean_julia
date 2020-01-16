README.txt for loading and plotting NPP CTD data

Revision: 2019-12-12 by Donglai Gong

NOTE: please follow the installation and setup steps exactly.

Purpose:
This set of Julia scripts will first loads CTD data in CNV format from the NPP expedition for additional data analysis, then plot against WOA data for profile comparison.

Installation:
0) Grab all NPP CTD data in CNV format and place into a data directory (e.g. /Users/gong/Research/NPP/data/cnv or c:\Users\rivera\Research\NPP\data\cnv)
1) Unzip all content in NPP.zip into its own directory (e.g. /Users/gong/Research/NPP/ or c:\Users\rivera\Research\NPP\), this is your work directory.
#2) Install Anaconda3. Open a new Terminal window, then type "conda update conda", and "conda update anaconda".
3a) Type "pip install cmocean" in the Terminal.
3b) Type "conda install scikit-learn".
4) Download and install Julia from https://julialang.org/ on to your computer. Version 1.3 or higher.

Setup:
1) Launch Julia.
2) At the julia prompt (julia>) type cd(ENV["HOME"] * "/Research/NPP/"), which brings you into your work directory.
3) Once your are in your work directory (type pwd() to verify), type include("npp_pkg_setup.jl"). this will install the necessary toolboxes for your Julia installation. You should only need to run this once, but it'll take a while to complete.
4) [Do this only if you have installed Anaconda separate on a Mac of Linux machine!!!] Move the included startup.jl file into the "/Users/Shabangin/.julia/config/", using Mirella's computer as example. You may need to create this directory using the command "mkdir -p /Users/Shabangin/.julia/config", then "mv startup.jl /Users/Shabangin/.julia/config/" to move the starup.jl file. Please note: That '.' in front of 'julia' is crucial.


Usage:
0) Check and modify if necessary the file "npp_path_setup.jl" in your work directory if needed so that its variables are consistent with how you set up your work directory and where you placed your data, as you have done above during installation. should only need to do this once.
1) Open the file "npp_flags.jl", make sure that both 'reloadflag' and 'saveflag' are set to 1 before you run "run_npp_load_ctd.jl" for the first time. for subsequent runs, once the JLD2 file is created, you can change 'reloadflag' to 0 if you want to. it's optional.
2a) Run the script "run_npp_load_ctd.jl" by typing include("run_npp_load_ctd.jl") from the julia> prompt.
2b) If everything ran as expected, then you should have three variables ('npp', 'ctd', 'transect') containing the CTD data loaded into julia now.
2c) For example: To access the depth data from the the first cast, type 'ctd[1].z'. to access the transect names, type 'npp.sectname'.
3) If you want to compare NPP CTD data with WOA data, then run "npp_plot_woa_ctd.jl". This will load the 'profile' variable that includes WOA profiles at NPP CTD stations. It'll also make the plots if the flag is set (in npp_flags.jl).
4) If you want to plot all the NPP CTD sections, run "npp_plot_sections.jl". Note, this might take a WHILE.
