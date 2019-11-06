# Background estimation steps

## Making the trees and setting up your area

Make and run the condor job on all MC and data:  

    ./running/runBatchJobs.py  -m TreeAnalyzer/background_estimation/predMacros/makeBETrees.C -b -i procdatasets.conf  -j jobs/4_30_beTreees -o out 

Compile the outputs into a single file per process:  

    ./running/compile.sh jobs/4_30_beTreees/ compiled/  

You will likely have to update the script to correspond to your actual datasets. Also, make the "compiled" directory.  

On your computer, run the following script to setup your area with the complete directory structure, download your files form lxplus, and skim the trees for use. You will have to update some hardcoded links in the script. This should be run in a clean directory that will be where you run all of the later scripts.  

    . background_estimation/predMacros/getAreaReady.sh /uscms/home/nmccoll/nobackup/2011-04-15-susyra2/rel_HbbWW/work/analyzer_running/compiled/ trees/ background_estimation/predMacros/skimTree.C  

__Parameters__

1. Directory where the compiled trees are
2. Where the trees will go. By default this should be "trees/"
3. Skimming script, the location you give the script should be accessible 

## Running the background estimation

This and all future steps assume that "background_estimation" is linked in your work directory (where you ran getAreaReady.sh).  

    . background_estimation/predMacros/runInputs.sh . background_estimation/predMacros/  

__Parameters__

1. Where to run
2. Where the macros are

It is that easy! Wait a few minutes and you will have all of your background and signal templates and you can make your datacards. It assumes that you have a few cores as there a lot of parallel processes running. Depending on what you want to do, you will probably want to turn some off. Here is the breakdown by block:  

### makeSignalInputs.C  
This macro makes all of the signal inputs.  
__Parameters__

1. Step of the estimation  
    * __0__ Make the histograms to be fit and the yields
    * __1__ Fit the histograms and make templates
    * __4__ Empty, just compile
2. Signal type enum
3. Where to find the trees
    
### makeTTBarSF.C
Make the ttbar normalization scale factor. There is a single scale factor for all regions  (control and search).  

__Parameters__    

1. Step of the estimation  
    * __0__ Make the histograms to be fit and then fit them
    * __1__ Test the scale factors
2. Where to find the trees
    
### makeBKGInputs.C
Make the background search region templates.

__Parameters__    

1. Enum of the background to run
2. Where to find the trees

For each background, the macro goes step by step to make the templates. In general, the steps are:

1. Preparation (like the QCD ratio)
2. Enum of the background to run
3. Where to find the trees

### Notes  

The last entry is to make the data distributions. (4) as the parameter is used for pseuod-data, while (-1) is used for the real data. At every step of the way there are a bunch of checks that can be run to verify that everything is working.

## Testing the templates  
Every step of the template production can be tested. This is done across a few different files. Each take three command line arguments:  

1. Step of the test (the tests go in the order of the makeBKGInputs, but the numbers are not 1-to-1 with those numbers)
2. Background or signal to test
3. For the background tests only: which region to test (0: SR, 1: tt CR, 2: q/g CR)
4. (3 for signal tests) directory name for the outputs  

The last one is important. While the the scripts will display the plots on your screen, if you give the script a name of a directory (that already exists), the pdfs and the root files of the tests will be stored. 

    rr 'background_estimation/plotting/plotSignalTests.C+(4,0,false ,"baseline")'  
    rr 'background_estimation/plotting/plotNonResBkgTests.C+(0,true,0,"baseline")'  
    rr 'background_estimation/plotting/plotResBkgTests.C+(0,true,0,"baseline")'   

### plotSignalTests.C  
Tests the signal, one thing to point out is that the signal yield tests also output the yields saved in the supplementary material.  
### plotNonResBkgTests.C  
Tests the q/g and lost t/w backgrounds.  
### plotResBkgTests.C  
Tests the mt and mW backgrounds. Step 6 is used to test all four backgrounds combined.  

## Make data cards  
Make the cards by running the following in the same directory as where you did the template building:  

    . background_estimation/predMacros/makeCards.sh nonCond_ePTRatio0p4_wQCD_sepQGemu 0 0 background_estimation/predMacros/

__Parameters__  

1. Card label, should be something unique so you can keep track of things.  
2. What region to run in (0 SR, 1 tt CR, 2 q/g CR).  
3. What type of signal you are running.  
4. Where the macros are located.  

The macro will also automaticly send the cards to fermilab, you may want to change where they do that or turn that off. This macro calls `makeCard.C`. This file is where the backgrounds and systematics are defined. Also, this file is where you say if you want to look at real or pseudo data.  

After running this you should have some data cards in `limits/CARDLABEL`. In particular, the `combinedCard.txt` and `combined.root` are the files corresponding to the card for all combined regions. The directory `plots` in here is where all of the plots go when you test these cards.  

## Test the model with data (and extra plots)  
Most of the tests are done with `plotDataTests.C`. We pretty much assume that you have a complete card at this point, but it will also use information from the inputs directories.  

### plotDataTests.C  
    rr `background_estimation/plotting/plotDataTests.C+(0,0,"limits/c_baseline2_data")' 

__Parameters__ 

1. Step of the test  
    *__0__ Run the post-fit so you can make post-fit plots (takes some time)
    *__1__ Test the pre-fit model
    *__2__ Test the post-fit model, as used in the AN
    *__3__ Do the saturated GOF test for projections of the model
    *__4__ Make summary plots, this is the global GOF test and the systematic pulls. This requires that you have run GoodnessOfFit and FitDiagnostic steps of the statistical tests.
    *__5__ Plot the bias test results. Requires that you run the bias tests.
    *__6__ Post-fit CR plots for the paper
    *__7__ Post-fit SR plots for the paper
    *__8__ Supplemental post-fit SR plots
    *__9__ More Supplemental post-fit SR plots
2. What region you are running
3. Data card directory

### plotDataTests.C  
This macro makes the search region variable plots for the paper. You can turn off blinding in the search region if you want. It runs in two steps, make the histograms in the first and the plots in the second. The second argument is for which region you are running on. The 1 and 2 steps have to be done after the 0 completes.

    rr -b -q 'background_estimation/plotting/plotSRVariables.C+   (0,0,"trees/betrees_mc.root","mc")' &  
    rr -b -q 'background_estimation/plotting/plotSRVariables.C+   (0,0,"trees/betrees_data.root","data")' &  
    rr -b -q 'background_estimation/plotting/plotSRVariables.C+   (0,0,"trees/out_radion_hh_bbinc_m1000_0.root","m1000")' &  
    rr -b -q 'background_estimation/plotting/plotSRVariables.C+   (0,0,"trees/out_radion_hh_bbinc_m2500_0.root","m2500")' &  
    rr 'background_estimation/plotting/plotSRVariables.C+         (1,0,"","")'   
    rr 'background_estimation/plotting/plotSRVariables.C+         (2,0,"","")'   

## Running limits and statistical tests  

### Asymptotic limits  
    nohup combine -m 800   -M AsymptoticLimits --run expected  --rMax 0.5 -v 2 combinedCard.txt &  
__Parameters__  

* __`nohup` and `&`__ Good to run in the background
* __`-m`__ Mass of the signal you want to test
* __`--run expected`__ Only run the expected limits...turn this off if you want the observed too
* __`--rMax 0.5`__ Maximum signal strength value...this should be tuned for the tested mass
* __`-v 2`__ Decent verbosity
* __`combinedCard.txt`__ your datacard 

The macro `limit_plotting/doLXPLimits.py` is used to do the limits as batch jobs at LXPlus.  

### Asymptotic observation significance

    combine -M Significance  -m 2300 ./combined.root

### Observation significance with toys  
This takes some time, so it we use crab and the combine tool to do it:  

    combineTool.py -d combinedCard.root -m 2500  -M HybridNew --LHCmode LHC-limits -v 2 --singlePoint 0.002:0.05:0.001 -T 100 --clsAcc 0 --iterations 2 --fork 2    --saveToys --saveHybridResult  --job-mode crab3 --task-name grid-test  --custom-crab cuscrab.py  
You care about `--singlePoint 0.002:0.05:0.001 -T 100 --clsAcc 0 --iterations 2 --fork 2` which is how the job is defined. The first argument is the scan of signal strength (min,max,binning). `-T 100` says do 100 toys. `--iterations 2` says do two batches of 100 toys. `--fork 2` says fork this into two process, you need some forking to get root to cleanup after each iteration. The problem is that there are memory leaks after each toy, so you have to balance between startup costs and the memory leak. You need `--custom-crab cuscrab.py` to tell the tool how to run on crab. The file should be:  

__cuscrab.py contents:__   

    def custom_crab(config):
        config.Site.storageSite = 'T3_US_FNALLPC'  

Then, after everything runs, combine all of the outputs:  

    cd crab_grid-test/results/
    for a in `ls -1 *.tar`; do tar -xvf $a; done
    cd ../.. 
    hadd higgsCombine.Test.HybridNew.m3500.123456.root crab_grid-test/results/higgsCombine*root  

Then get the final result:  

    combine -d combinedCard.root -m 2500  -M HybridNew --LHCmode LHC-limits -v 2 --readHybridResults --grid=higgsCombine.Test.HybridNew.mH2500.123456.root  --expectedFromGrid=0.5  

You can also do the same thing....but interactively. This time you need to break up your iterations and give it random seeds. An example:  

    combine -M HybridNew -m 2300 ../combined.root --LHCmode LHC-significance  --saveToys --fullBToys --saveHybridResult -T 10 -i 5 -s -1 &
    combine -M HybridNew -m 2300 ../combined.root --LHCmode LHC-significance  --saveToys --fullBToys --saveHybridResult -T 10 -i 5 -s -1 &
    combine -M HybridNew -m 2300 ../combined.root --LHCmode LHC-significance  --saveToys --fullBToys --saveHybridResult -T 10 -i 5 -s -1 &
    combine -M HybridNew -m 2300 ../combined.root --LHCmode LHC-significance  --saveToys --fullBToys --saveHybridResult -T 10 -i 5 -s -1 &
    combine -M HybridNew -m 2300 ../combined.root --LHCmode LHC-significance  --saveToys --fullBToys --saveHybridResult -T 10 -i 5 -s -1 &
    combine -M HybridNew -m 2300 ../combined.root --LHCmode LHC-significance  --saveToys --fullBToys --saveHybridResult -T 10 -i 5 -s -1 &
    combine -M HybridNew -m 2300 ../combined.root --LHCmode LHC-significance  --saveToys --fullBToys --saveHybridResult -T 10 -i 5 -s -1 &
    hadd higgsCombineTest.HybridNew.mH2300.allToys.root higgsCombineTest.HybridNew.mH2300.*    
    combine -M HybridNew -m 2300 ../combined.root --LHCmode LHC-significance   --readHybridResult --toysFile=higgsCombineTest.HybridNew.mH2300.allToys.root 

### Do a fit and comput the standard fit diagnostics  
These are a part of the "standard tests." You will need to make sure that you have your card in root file form:  

    text2workspace.py combinedCard.txt -o combined.root

Do a fit to the data with some mass hypothesis:  

    combine -M MaxLikelihoodFit -m 1000 combined.root  

Get the fit diagnostics (e.g. NP pulls):  

     python ~/GitRepositories/TreeAnalyzer/framework/HiggsAnalysis/CombinedLimit/test/diffNuisances.py -g plots.root fitDiagnostics.root 

### Goodness of fit tests  
Run first, to get the observed value:

    combine -M GoodnessOfFit combinedCard.txt --algo saturated --fixedSignalStrength 0  -m 1000 --toysFreq

Then run your toys:

    combine -M GoodnessOfFit combinedCard.txt --algo saturated --fixedSignalStrength 0  -m 1000 --toysFreq -t 100 
    combine -M GoodnessOfFit combinedCard.txt --algo saturated --fixedSignalStrength 0  -m 1000 --toysFreq -t 100 -s 1 &
    combine -M GoodnessOfFit combinedCard.txt --algo saturated --fixedSignalStrength 0  -m 1000 --toysFreq -t 100 -s 2 &
    combine -M GoodnessOfFit combinedCard.txt --algo saturated --fixedSignalStrength 0  -m 1000 --toysFreq -t 100 -s 3 &
    combine -M GoodnessOfFit combinedCard.txt --algo saturated --fixedSignalStrength 0  -m 1000 --toysFreq -t 100 -s 4 &
    combine -M GoodnessOfFit combinedCard.txt --algo saturated --fixedSignalStrength 0  -m 1000 --toysFreq -t 100 -s 5 &
    combine -M GoodnessOfFit combinedCard.txt --algo saturated --fixedSignalStrength 0  -m 1000 --toysFreq -t 100 -s 6 &
    combine -M GoodnessOfFit combinedCard.txt --algo saturated --fixedSignalStrength 0  -m 1000 --toysFreq -t 100 -s 7 &
    combine -M GoodnessOfFit combinedCard.txt --algo saturated --fixedSignalStrength 0  -m 1000 --toysFreq -t 100 -s 8 &
    hadd  higgsCombineTest.GoodnessOfFit.mH2000.toys.root  higgsCombineTest.GoodnessOfFit.mH2000.*.root

### Calculate the systematic impacts  
Start by getting your input root file:  

    text2workspace.py combinedCard.txt -m 1000 -o forImpact.root

Step one is to do the initial fit, it depends on the mass and the test that you are doing.  

    combineTool.py -M Impacts -d forImpact.root  --doInitialFit --robustFit 1 -m 1000 --rMax 0.42 --rMin 0.047 --toysFrequentist -t -1 --expectSignal=0.14  

Do a 1 TeV signal with bounds on the fitted signal strength (important if you want it to converge). `--toysFrequentist -t -1` says use the Asimov dataset. `--expectSignal=0.14` say inject signal with a certain strength. You can change this to 0.  

    combineTool.py -M Impacts -d forImpact.root --doInitialFit --robustFit 1 -m 1000 --rMax 0.2  

Is very similar, byt this time you are fitting to data. Then you have to do the second step...actually computing the impacts. How you do the fit and the toys have to be the same as step one. If you wanted to run interactively for the previous two tests you would:  

    combineTool.py -M Impacts -d forImpact.root --doFits --robustFit 1 --parallel 20 -m 1000 --rMax 0.42 --rMin 0.047 --toysFrequentist -t -1 --expectSignal=0.14
    combineTool.py -M Impacts -d forImpact.root --doFits --robustFit 1 --parallel 20 -m 1000 --rMax 0.2

If you wanted to run with crab instead (crab.py is the same as for observation with toys):  

    combineTool.py -M Impacts -d forImpact.root --doFits --robustFit 1 --parallel 20 -m 1000 --rMax 0.42 --rMin 0.047 --toysFrequentist -t -1 --expectSignal=0.14 --job-mode crab3 --task-name grid-test  --custom-crab cuscrab.py 
    combineTool.py -M Impacts -d forImpact.root --doFits --robustFit 1 --parallel 20 -m 1000 --rMax 0.2 --job-mode crab3 --task-name grid-test  --custom-crab cuscrab.py 

Finally, you can then make your pretty plots:  

    plotImpacts.py -i impacts.json -o impacts

