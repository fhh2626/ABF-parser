# ABF-parser
A tool parsing and analyzing ABF outputs of NAMD

See /unitTest/unitTest.cpp for how to use this tool as a library. This tool can also be used as a standalone program. One can compile this tool by
```
g++ main.cpp
```
Then one can simply run
```
main.exe config.ini 0
```
The third parameter, 0-6, defines the job type:
```
enum JobType {
        kCalcError,
        kFindPathway,
        kCalcRMSD,
        kCalcFreeEnergyDiff,
        kCalcEnergyInCVSpace,
        kMergeWindows,
        kReweightPmf,
        kProjectPmf
    };
```
The second parameter indicates a config file. See /unitTest for examples about how to write a config file.
