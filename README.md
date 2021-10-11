# MiniAODBphysicsForRun3

## Intructions to run the examples.
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_12_1_0_pre2
cd CMSSW_12_1_0_pre2/src/
cmsenv
voms-proxy-init -voms cms -valid 192:00
git clone https://github.com/jmejiagu/MiniAODBphysicsUltraLegacyRun2.git myAnalyzers/JPsiKsPAT
cd myAnalyzers/JPsiKsPAT/
git checkout master
cd ../..
scram b

```

Run: (use your favorite input sample. You will see examples in the confi files)


```
cd myAnalyzers/JPsiKsPAT/test/
cmsRun PsikaonRootupler.py
cmsRun Psikaon_MC_Rootupler
```

This example is for Bu hadron, both data and MC config files. In test directory you can find other examples to run Bs and Bd hadrons too.