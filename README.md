# ChicAnalysisScripts
Collection of python scripts for analysing chic 

## Executing macros 
To run the python macros I used the following steps:
1) Give the function in the macro the correct input and save the file
2) enter O2Physics environment
3) go to the folder where the macro is located
4) run python name_macro.py

## General structure of the macros
In all macros the structure of the function for the histograms is the same. The input to the function consists of the input root file, a number for the option and the name of the file where the histogram is saved. First the data necessary for the histograms are collected from the input root file. Then a TCanvas window is created and the area of the histogram is defined. Optional a calculation with the data from the root file takes place. The Histograms x- and y-axis labels are set. In if clauses the histograms are filled with the data and the legend is specified. Afterwards the histogram gets some additional text (Simulation this thesis, pp, sqrt(s)=13.6TeV). The last step is the saving of the histogram in the given input string. There are some comments in the codes to what the code lines do. 

## Details to the macros
The **acceptance.py** file has the cutacceptance_plot() function. This function calculates the acceptance for J/psi -> e^{+}e^{-} (option 0), e^{+}e^{-} from chic1 and chic2 (option 12) and chic1 and chic2 -> photon+e^{+}e^{-} (option 7). The acceptance histograms have variable binning: [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.,2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22].

The **varefficiency.file** has the cutefficiency_plot function which determines the efficiency with variable binning by dividing the MC matched candidates by the MC generated true candidates in acceptance.
- option 0: J/psi -> e^{+}e^{-}  
- option 11: e^{+}e^{-} from chic1 and chic2 
- option 20: photon from chic1 and chic2
- option 17: chic1 and chic2 -> photon e^{+} e^{-} 
- binning for option 0 and 11: [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.,2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20]
- binning for option 20: [0, 0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.5,3,3.5,4.,4.5]
- binning for option 17: [0, 0.5, 1, 1.5, 2., 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20]

The **deltaMass.py** file has the deltamass_plot() function to show the delta mass histograms with arrows at the locations of chic1 and chic2 PDG masses.
- option 13: deltamass of photon e^{+} e^{-} MC reconstructed
- option 23: deltamass of photon e^{+} e^{-} MC generated true
- option 26: deltamass of photon e^{+} e^{-} MC reconstructed matched
- option 33: deltamass of photon e^{+} e^{-} MC reconstructed and MC reconstructed matched
- option 50: deltamass of photon e^{+} e^{-} and mass of photon e^{+} e^{-} MC reconstructed matched
There are more options which can be used to draw the delta mass histograms of other combinations.

The **itdeltaMass.py** file has the deltamass_plot() function which can fit a Breit-Wigner or an asymmetric Gaussian to the delta mass distributions of MC generated true and MC reconstructed matched triple candidates.
- option 0: Breit-Wigner fit to MC generated true triple candidates
- option 4: asymmetric Gaussian fit to MC reconstructed matched triple candidates
For new data the parameter limits and parameters input values have to probably be adjusted.

The **etapT.py** file has the function pT_plot(). It draws all the p_T functions. 
- option 0: chic1 and chic2 MC generated true
- option 3: J/psi and J/psi from chic1 and chic2 MC generated true
- option 6: J/psi -> ee and ee from chic1 and chic2 MC generated true
- option 13: J/psi -> ee and ee from chic1 and chic2 MC reconstructed matched
- option 18: chic1 and chic2 -> photon e^{+} e^{-} MC reconstructed matched
- option 27: photon from chic1 and chic2 MC reconstructed matched and MC generated true

The **mass.py** has the function mass_plot() which draws all the mass histrograms with arrows at PDG J/psi, chic1 and chic2 mass locations.
- option 0: J/psi -> e^{+} e^{-} MC reconstructed matched
- option 3: chic1 and chic2 -> photon e^{+} e^{-} MC reconstructed
- option 4: chic1 and chic2 -> photon e^{+} e^{-} MC reconstructed and MC reconstructed matched
- option 5: chic1 and chic2 -> photon e^{+} e^{-} MC reconstruced matched
- option 10: J/psi -> e^{+} e^{-} MC generated true
- option 20: e^{+} e^{-} and e^{+} e^{+} + e^{-} e^{-} MC reconstructed
- option 21: e^{+} e^{-} minus e^{+} e^{+} + e^{-} e^{-} MC reconstructed and J/psi -> e^{+} e^{-} MC reconstructed matched
- option 22: e^{+} e^{-} MC reconstructed small range
- option 30: chic1 and chic2 -> photon e^{+} e^{-} MC generated true

The **photons.py** file has the photon_plot function which creates a pT distribution plot, 2D photon conversion points and photon conversion probability calculation: 
- option 0: pT distribution all photons, converted photons and reconstructed photons
- option 1: conversion probability
- option 4: photon conversion points MC generated true
- option 5: photon conversion points MC reconstructed
- option 6: photon conversion points MC reconstructed matched to true conversion coordinates

The **ratiophotons.py** file has the ratiophoton_plot() function which calculates the number of photon conversion per conversion radius and efficiency of reconstructing photons depending on the conversion radius.
- option 1: number of photon conversion per conversion radius for MC generated true, MC reconstructed and MC reconstructed matched to true conversion coordinates
- option 2: efficiency of reconstructing photons depending on conversion radius by dividing number of photon conversions MC reconstructed matched to true conversion coordinates by MC generated true

The **cuthistograns.py** file has the cut_plot() function which plots histogram containing the cut parameters.
- option 1: cosine of pointing angle  
- option 2: dE/dx signal of electrons from PCM
- option 3: Armenteros-Podolanski plot consisting only photons
- option 4: dE/dx signal of primary photons

## Configuration file for AnalysisDileptonPhoton
The **configuration file (config_AnalysisDileptonPhoton_dqEfficiency.json)** is for running the AnalysisDileptonPhoton in dqEfficiency with MC data.

This is the command for running over MC data:
o2-analysis-dq-table-maker-mc --configuration json://tempConfigMCchic.json --severity error --shm-segment-size 12000000000 --aod-writer-json aodWriterTempConfig.json -b | 
o2-analysis-timestamp --configuration json://tempConfigMCchic.json -b | 
o2-analysis-event-selection --configuration json://tempConfigMCchic.json -b | 
o2-analysis-multiplicity-table --configuration json://tempConfigMCchic.json -b | 
o2-analysis-trackselection --configuration json://tempConfigMCchic.json -b | 
o2-analysis-pid-tof-base --configuration json://tempConfigMCchic.json -b | 
o2-analysis-pid-tof --configuration json://tempConfigMCchic.json -b | 
o2-analysis-pid-tof-full --configuration json://tempConfigMCchic.json -b | 
o2-analysis-pid-tof-beta --configuration json://tempConfigMCchic.json -b | 
o2-analysis-pid-tpc-base --configuration json://tempConfigMCchic.json -b | 
o2-analysis-pid-tpc-full --configuration json://tempConfigMCchic.json -b | 
o2-analysis-track-propagation --configuration json://tempConfigMCchic.json -b  | 
o2-analysis-dq-efficiency --configuration json://tempConfigMCchic.json -b |
o2-analysis-em-photon-conversion-builder -b --configuration json://tempConfigMCchic.json |
o2-analysis-em-create-emreduced-event -b --configuration json://tempConfigMCchic.json | 
o2-analysis-em-associate-mc-info -b --configuration json://tempConfigMCchic.json | 
o2-analysis-em-pcm-qc-mc -b --configuration json://tempConfigMCchic.json 

The **configuration file (config_AnalysisDileptonPhoton_tableReader.json)** is for running AnalysisDileptonPhoton in the tableReader.cxx for ALICE Run 3 data. 

This is the command for running over ALICE Run 3 data: 
o2-analysis-pid-tof-full -b --configuration json://configuration-DQ.json | 
o2-analysis-pid-tof-beta -b --configuration json://configuration-DQ.json | 
o2-analysis-pid-tof-base -b --configuration json://configuration-DQ.json | 
o2-analysis-pid-tof -b --configuration json://configuration-DQ.json | 
o2-analysis-ft0-corrected-table -b --configuration json://configuration-DQ.json | 
o2-analysis-timestamp -b --configuration json://configuration-DQ.json | 
o2-analysis-pid-tpc-full -b --configuration json://configuration-DQ.json | 
o2-analysis-pid-tpc-base -b --configuration json://configuration-DQ.json | 
o2-analysis-track-propagation -b --configuration json://configuration-DQ.json | 
o2-analysis-multiplicity-table -b --configuration json://configuration-DQ.json | 
o2-analysis-trackselection -b --configuration json://configuration-DQ.json | 
o2-analysis-tracks-extra-converter -b --configuration json://configuration-DQ.json | 
o2-analysis-event-selection -b --configuration json://configuration-DQ.json | 
o2-analysis-dq-table-maker -b --configuration json://configuration-DQ.json | 
o2-analysis-dq-table-reader -b --configuration json://configuration-DQ.json | 
o2-analysis-em-create-pcm -b --configuration json://configuration-DQ.json | 
o2-analysis-em-create-emreduced-event -b --configuration json://configuration-DQ.json | 
o2-analysis-em-skimmer-gamma-conversion -b --configuration json://configuration-DQ.json | 
o2-analysis-em-pcm-qc -b --configuration json://configuration-DQ.json --aod-file 

