I use triggerCondor scripts to run Fun4All, create the root files, and then use makeHistogram scripts to make histograms for each root file. 
And then combine hisograms to combine them for a single run.

makeFileList takes all the files in a specified directory and turns into a list file, I use this to feed into condor submit file for Histogram generation after Fun4All step
