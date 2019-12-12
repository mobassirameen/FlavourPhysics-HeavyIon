# BsToJpsiPhi
The BsToJpsiPhi is the reconstruction code for the Run-II Bs->Jpsiphi( Jpsi->mumu, Phi->KK) analysis, 
works to collect the muons and packed particle flow objects with proper trigger level object selection. The gloal is to 
collect all kinematic observables for the flavour tagged agular analysis. The kinematic observables are three angles in the 
transversity basis of the decay process and the lifetime and lifetime uncertainty of the Bs mother candidate.

#The BsToJpsiPhi/src directory 

This is main C++ analyzer module which is further wrapped with a python script resides in the test directory along 
with some small sripts ( viz., crab script)

# Codes 

There are several C++ macros which are essentially used for the different studies such as efficieny calculation, model the 
mass and ct-uncertainty, model 2d mass-mistag using moment morphing procedure. 
