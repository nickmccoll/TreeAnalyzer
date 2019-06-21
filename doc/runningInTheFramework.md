##Running a macro

All macros take the same command line arguments, this allows them to be run at the command line or in the batch system:  

    rr -b -q 'background_estimation/predMacros/makeBETrees.C+("FILE",1,1,"test.root")'  
__Parameters:__  

1.  Input filename
2.  Tree integer (0 for data, 1 for MC)
3.  Random number generator seed
4.  Output filename
5.  Sample cross section (does not need to be filled)
6.  Total number of events in the sample (does not need to be filled)

The cross section and number of events is used to normalize the MC samples. This is usually stored in a configuration file.

Example file:  
/eos/uscms/store/user/lpchh/HHWWbb_trees/3_21_19_complete/BulkGravTohhTohVVhbb_narrow_M-1200_TuneCP5_13TeV-madgraph-pythia8/crab_BulkGravTohhTohVVhbb_narrow_M-1200_TuneCP5_13TeV-madgraph-pythia8/190322_041954/0000/BulkGravTohhTohVVhbb_narrow_M-1200_TuneCP5_13TeV-madgraph-pythia8_1.root



##Running from a configuration file
    ./running/runBatchJobs.py  -m TreeAnalyzer/background_estimation/predMacros/makeBETrees.C -b -i procdatasets.conf  -j jobs/4_30_beTreees -o out 

__Parameters:__  

* -m : The macro to run
* -b : Use the batch system (remove to run interactively)
* -i : configuration file telling you which samples to run on
* -j : job directory to make all of the temporary files needed to setup the jobs
* -o : output directory to store job output. If you run interactively the output goes here, if not it goes to the job dir

##Configuration files

Configuration files are stored under the running directory (eg. procdatasets_2017.conf)

    TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8			1	364.3508	99469578	20	/eos/uscms/store/user/lpchh/HHWWbb_trees/2_12_19_complete/

__Column descriptions__

1. Sample name used to find the input files in the input dir and to label the output files
2. Tree integer (0 for data, 1 for MC)
3. Sample cross section (in pb)
4. Total number of events in the sample
5. Number of files to process per job
6. Directory to find these files. The script will recursively look for all root files with the sample name.

