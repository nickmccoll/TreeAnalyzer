# TreeAnalyzer

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