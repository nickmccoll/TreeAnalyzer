## Setup commands
Execute the following commands on your system. If you are running at FNAL you will need root...set this up with a `cmsenv`.


```Shell
git clone git@github.com:{USERNAME}/TreeAnalyzer.git
cd TreeAnalyzer
git submodule init
git submodule update
cd framework
make
```

### Optional steps for combine-dependent macros
If you want to use some background estimation macros that depend on combine, you have to compile and setup combine. In the case of using linux (with CMSSW) you can use the built in combine...otherwise you have to install it.

#### CMSSW based instructions
```Shell
cd $CMSSW_BASE/src/
cmsenv
cd {TREEANALYZER_PATH}/framework/HiggsAnalysis/CombinedLimit
make
```

#### MacOS based instructions
Download and unzip Boost (https://www.boost.org) to some directory `{BOOSTSRC}`. Then follow the following installation instructions in some build directory `{BOOSTBUILD}`:

```Shell
mkdir {BOOSTBUILD}
cd {BOOSTSRC}
./bootstrap.sh --prefix={BOOSTBUILD}
./b2 install
```
Boost does not have good inter-library linking for MacOS so you will have to additionally do the following:

```Shell
cd {BOOSTBUILD}/lib
vim change.sh
---- COPY AND PASTE THE SCRIPT IN THE NEXT BLOCK ----
chmod +x change.sh
./change.sh
```

Where you copy and paste the following script:

```Shell
#!/bin/bash
  
# Modify the absolute dylib paths baked into the libraries
for i in *.dylib
do
    FULLPATH=`pwd`/$i
    install_name_tool -id $FULLPATH $i
    echo -change $i $FULLPATH
    done > changes
for i in *.dylib
do
    install_name_tool `cat changes` $i
done
rm changes
```

Now change your bash profile by adding the following lines:
```Shell
export BOOSTPATH="{BOOSTBUILD}"
export HIGGSCOMBPATH="{TREEANALYZER_PATH}/framework/HiggsAnalysis/CombinedLimit/"
export PYTHONPATH=${PYTHONPATH}:${HIGGSCOMBPATH}/lib/python:${HIGGSCOMBPATH}/lib
export PATH=${PATH}:${HIGGSCOMBPATH}/scripts
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:${HIGGSCOMBPATH}/lib:${BOOSTPATH}/lib
```

Where `{BOOSTBUILD}` is to be replaced with what you used above and `{TREEANALYZER_PATH}` is path that you downloaded the repo. Make sure that you open a new terminal window so that you have a good enviornment and install combine:

```Shell
cd {TREEANALYZER_PATH}/framework/HiggsAnalysis/CombinedLimit
make
```

## Enviornment
Create the following `rootLogon.C` file.  Here `{TREEANALYZER_PATH}` is path that you downloaded the repo. For example,  `/uscms/home/nmccoll/nobackup/2011-04-15-susyra2/rel_HbbWW/TreeAnalyzer`. Is where I put mine.

```C++
{
    gSystem->AddIncludePath(" -I{TREEANALYZER_PATH}/framework ");
    gSystem->Load("{TREEANALYZER_PATH}/framework/libTreeAnalysis.so");
}
```  

You will also need a `.rootrc` file telling root to load this `rootLogon.C`.  For example, you can place this in your working directory:  

```Shell
Rint.Logon:              rootlogon.C
Rint.Logoff:             rootlogoff.C
Unix.*.Root.MacroPath:    .:$(HOME):
```

### Optional steps for combine-dependent macros
If you want to use the macros that depend on combine add the following to your `rootLogon.C` file:

```Shell
gSystem->Load("{TREEANALYZER_PATH}/framework/HiggsAnalysis/CombinedLimit/lib/libHiggsAnalysisCombinedLimit");
```

You will also need to make a `.rootlogon.py` file in your home directory. Fill it with:

```Shell
import ROOT
ROOT.gSystem.AddDynamicPath("{TREEANALYZER_PATH}/framework/HiggsAnalysis/CombinedLimit/lib")
```

## Data directory
Download the data directory:  

```Shell
mkdir data
cd data
wget -O data.zip  https://www.dropbox.com/sh/xo9djqy4trdtoww/AAAKKZK4xez-ICT09xgIDu7Aa?dl=1
unzip data.zip
rm data.zip
cd ..
```

You can download it into your work directory, jobs will run automatically if they are run with the "data" directory in the execution directory. Otherwise you can set an enviornment variable to pick it up from some other location: `export TREEANALYZER_DATA="/Users/nmccoll/Dropbox/Work/Projects/HHbbWW/hbbww_data/"`
