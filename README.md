# MiniAODBphysicsForRun3

## Intructions to run the examples.
```
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc900
cmsrel CMSSW_12_1_0_pre2
cd CMSSW_12_1_0_pre2/src/
cmsenv
voms-proxy-init -voms cms -valid 192:00
git clone git@github.com:jmejiagu/SoftandFrame_MiniAOD.git myAnalyzers/JPsiKsPAT
cd myAnalyzers/JPsiKsPAT/
git checkout main
cd ../..
scram b

```

Run: (use your favorite input sample. You will see examples in the confi files)


```
cd myAnalyzers/JPsiKsPAT/test/
cmsRun PsikaonRootupler.py
cmsRun Psikaon_MC_Rootupler
```

This example is for JpsiTrack reconstruction. In particular, the config file refers to Bu hadron, both data and MC config files. 
In test directory you can find examples to run JpsiTrkTrk, JPsiV0, and others.

Run Chi(c,b)

```
cd myAnalyzers/JPsiKsPAT/test/
cmsRun run-chic-miniaod.py
cmsRun run-chib-miniaod.py
```

