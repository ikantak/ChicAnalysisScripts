import re
import numpy as np
import datetime
import math

import ROOT
import ctypes
from ctypes import *
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency, TArrow
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX

gStyle.SetOptStat(0);
gStyle.SetOptTitle(0);



def make_common_style(g1,marker,size,color,width=1,fill=0):
    g1.SetMarkerStyle(marker);
    g1.SetMarkerColor(color);
    g1.SetMarkerSize(size);
    g1.SetLineColor(color);
    g1.SetLineWidth(width);
    g1.SetFillColor(color);
    g1.SetFillStyle(fill);

def deltamass_plot(filename1,  option, plotname):
    rootfile_data = TFile.Open(filename1, "READ");
    list_data = rootfile_data.Get("analysis-dilepton-photon");
    list_data2 = list_data.Get("output");
     
    if (option == 0 or option == 16): 
        mc_eephoton_chic1 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic1")
        mc_eephoton_chic2 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic2")
        eephoton1_deltaM = mc_eephoton_chic1.FindObject("DeltaMass");
        eephoton1_deltaM.SetDirectory(0);
        make_common_style(eephoton1_deltaM, kFullCrossX, 1.0, kSpring-8, 1, 0);
        ROOT.SetOwnership(eephoton1_deltaM, False);
        eephoton1_deltaM.Scale(1, "width")

        eephoton2_deltaM = mc_eephoton_chic2.FindObject("DeltaMass");
        eephoton2_deltaM.SetDirectory(0);
        make_common_style(eephoton2_deltaM, kFullCrossX, 1.0, kViolet+2, 1, 0);
        ROOT.SetOwnership(eephoton2_deltaM, False);
        eephoton2_deltaM.Scale(1, "width")
    
    if option == 1: 
        mc_eephoton_chic12 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic12")
        eephoton12_deltaM = mc_eephoton_chic12.FindObject("DeltaMass");
        eephoton12_deltaM.SetDirectory(0);
        make_common_style(eephoton12_deltaM, kFullCrossX, 1.0, kTeal-1, 1, 0);
        ROOT.SetOwnership(eephoton12_deltaM, False);
        eephoton12_deltaM.Scale(1, "width")

    if (option == 2 or option == 20): 
        mc_eephoton_chic1 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic1")
        mc_eephoton_chic2 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic2")
        mc_eephoton_chic12 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic12")
        eephoton1_deltaM = mc_eephoton_chic1.FindObject("DeltaMass");
        eephoton1_deltaM.SetDirectory(0);
        make_common_style(eephoton1_deltaM, kFullCrossX, 1.0, kRed, 1, 0);
        ROOT.SetOwnership(eephoton1_deltaM, False);
        eephoton1_deltaM.Scale(1, "width")

        eephoton2_deltaM = mc_eephoton_chic2.FindObject("DeltaMass");
        eephoton2_deltaM.SetDirectory(0);
        make_common_style(eephoton2_deltaM, kFullCrossX, 1.0, kOrange+7, 1, 0);
        ROOT.SetOwnership(eephoton2_deltaM, False);
        eephoton2_deltaM.Scale(1, "width")
    
        eephoton12_deltaM = mc_eephoton_chic12.FindObject("DeltaMass");
        eephoton12_deltaM.SetDirectory(0);
        make_common_style(eephoton12_deltaM, kFullCrossX, 1.0, kRed+2, 1, 0);
        ROOT.SetOwnership(eephoton12_deltaM, False);
        eephoton12_deltaM.Scale(1, "width")
    
    if (option == 3 or option == 30): 
        data_dileptonphoton = list_data2.FindObject("DileptonPhotonInvMass_cut")
        dileptonphoton_deltaM = data_dileptonphoton.FindObject("DeltaMass");
        dileptonphoton_deltaM.SetDirectory(0);
        make_common_style(dileptonphoton_deltaM, kFullCrossX, 1.0, kMagenta+2, 1, 0);
        ROOT.SetOwnership(dileptonphoton_deltaM, False);
        dileptonphoton_deltaM.Scale(1, "width")
    
    if option == 4:
        mc_eephoton_chic1 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic1")
        eephoton1_deltaM = mc_eephoton_chic1.FindObject("DeltaMass");
        eephoton1_deltaM.SetDirectory(0);
        make_common_style(eephoton1_deltaM, kFullCrossX, 1.0, kSpring-8, 1, 0);
        ROOT.SetOwnership(eephoton1_deltaM, False);
        eephoton1_deltaM.Scale(1, "width")

    if option == 5: 
        mc_eephoton_chic2 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic2")
        eephoton2_deltaM = mc_eephoton_chic2.FindObject("DeltaMass");
        eephoton2_deltaM.SetDirectory(0);
        make_common_style(eephoton2_deltaM, kFullCrossX, 1.0, kViolet+2, 1, 0);
        ROOT.SetOwnership(eephoton2_deltaM, False);
        eephoton2_deltaM.Scale(1, "width")
     
    
    if (option == 10 or option == 16 or option == 12 or option == 20): 
        data_eephoton_chic1 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic1")
        data_eephoton_chic2 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic2")
        eephoton1_data_deltaM = data_eephoton_chic1.FindObject("DeltaMass");
        eephoton1_data_deltaM.SetDirectory(0);
        make_common_style(eephoton1_data_deltaM, kFullCrossX, 1.0, kPink-2, 1, 0);
        ROOT.SetOwnership(eephoton1_data_deltaM, False);
        eephoton1_data_deltaM.Scale(1, "width")

        eephoton2_data_deltaM = data_eephoton_chic2.FindObject("DeltaMass");
        eephoton2_data_deltaM.SetDirectory(0);
        make_common_style(eephoton2_data_deltaM, kFullCrossX, 1.0, kPink+1, 1, 0);
        ROOT.SetOwnership(eephoton2_data_deltaM, False);
        eephoton2_data_deltaM.Scale(1, "width")
    
    if (option == 11 or option == 12 or option == 20 or option == 30): 
        data_eephoton_chic12 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic12")
        eephoton12_data_deltaM = data_eephoton_chic12.FindObject("DeltaMass");
        eephoton12_data_deltaM.SetDirectory(0);
        make_common_style(eephoton12_data_deltaM, kFullCrossX, 1.0, kPink-8, 1, 0);
        ROOT.SetOwnership(eephoton12_data_deltaM, False);
        eephoton12_data_deltaM.Scale(1, "width")
    
    if (option == 13 ):
        data_dileptonphoton = list_data2.FindObject("DileptonPhotonInvMass_cut")
        dileptonphoton_deltaM_jpsi = data_dileptonphoton.FindObject("DeltaMass_Jpsi");
        dileptonphoton_deltaM_jpsi.SetDirectory(0);
        make_common_style(dileptonphoton_deltaM_jpsi, kFullCrossX, 1.0, kPink+7, 1, 0);
        ROOT.SetOwnership(dileptonphoton_deltaM_jpsi, False);
        dileptonphoton_deltaM_jpsi.Rebin(5)
        print(dileptonphoton_deltaM_jpsi.GetNbinsX())
        dileptonphoton_deltaM_jpsi.Scale(1, "width")

    if option == 14: 
        data_eephoton_chic1 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic1")
        eephoton1_data_deltaM = data_eephoton_chic1.FindObject("DeltaMass");
        eephoton1_data_deltaM.SetDirectory(0);
        make_common_style(eephoton1_data_deltaM, kFullCrossX, 1.0, kSpring-8, 1, 0);
        ROOT.SetOwnership(eephoton1_data_deltaM, False);
        eephoton1_data_deltaM.Scale(1, "width")

    if option == 15: 
        data_eephoton_chic2 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic2")
        eephoton2_data_deltaM = data_eephoton_chic2.FindObject("DeltaMass");
        eephoton2_data_deltaM.SetDirectory(0);
        make_common_style(eephoton2_data_deltaM, kFullCrossX, 1.0, kViolet+2, 1, 0);
        ROOT.SetOwnership(eephoton2_data_deltaM, False);
        eephoton2_data_deltaM.Scale(1, "width")
    
    if (option == 21 or option == 23):
        mc_eephoton_chic1 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic1")
        mc_eephoton_chic2 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic2")
        eephoton1_deltaMjpsi = mc_eephoton_chic1.FindObject("DeltaMass_Jpsi");
        eephoton1_deltaMjpsi.SetDirectory(0);
        make_common_style(eephoton1_deltaMjpsi, kFullCrossX, 1.0, kRed, 1, 0);
        ROOT.SetOwnership(eephoton1_deltaMjpsi, False);
        eephoton1_deltaMjpsi.Scale(1, "width")

        eephoton2_deltaMjpsi = mc_eephoton_chic2.FindObject("DeltaMass_Jpsi");
        eephoton2_deltaMjpsi.SetDirectory(0);
        make_common_style(eephoton2_deltaMjpsi, kFullCrossX, 1.0, kOrange+7, 1, 0);
        ROOT.SetOwnership(eephoton2_deltaMjpsi, False);
        eephoton2_deltaMjpsi.Scale(1, "width")
    
    if (option == 22 or option == 23):
        mc_eephoton_chic12 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic12")
        eephoton12_deltaMjpsi = mc_eephoton_chic12.FindObject("DeltaMass_Jpsi");
        eephoton12_deltaMjpsi.SetDirectory(0);
        make_common_style(eephoton12_deltaMjpsi, kFullCrossX, 1.0, kRed+2, 1, 0);
        ROOT.SetOwnership(eephoton12_deltaMjpsi, False);     
        eephoton12_deltaMjpsi.Scale(1, "width")   
    

    if (option == 24 or option ==26): 
        data_eephoton_chic1 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic1")
        data_eephoton_chic2 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic2")
        eephoton1_data_deltaMjpsi = data_eephoton_chic1.FindObject("DeltaMass_Jpsi");
        eephoton1_data_deltaMjpsi.SetDirectory(0);
        make_common_style(eephoton1_data_deltaMjpsi, kFullCrossX, 1.0, kViolet-1, 1, 0);
        ROOT.SetOwnership(eephoton1_data_deltaMjpsi, False);
        eephoton1_data_deltaMjpsi.Rebin(5)
        print(eephoton1_data_deltaMjpsi.GetNbinsX())
        eephoton1_data_deltaMjpsi.Scale(1, "width")

        eephoton2_data_deltaMjpsi = data_eephoton_chic2.FindObject("DeltaMass_Jpsi");
        eephoton2_data_deltaMjpsi.SetDirectory(0);
        make_common_style(eephoton2_data_deltaMjpsi, kFullCrossX, 1.0, kViolet-4, 1, 0);
        ROOT.SetOwnership(eephoton2_data_deltaMjpsi, False);
        eephoton2_data_deltaMjpsi.Rebin(5)
        print(eephoton2_data_deltaMjpsi.GetNbinsX())
        eephoton2_data_deltaMjpsi.Scale(1, "width")
    
    if (option == 25 or option == 26 ): 
        data_eephoton_chic12 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic12")
        eephoton12_data_deltaMjpsi = data_eephoton_chic12.FindObject("DeltaMass_Jpsi");
        eephoton12_data_deltaMjpsi.SetDirectory(0);
        make_common_style(eephoton12_data_deltaMjpsi, kFullCrossX, 1.0, kMagenta+3, 1, 0);
        ROOT.SetOwnership(eephoton12_data_deltaMjpsi, False);
        eephoton12_data_deltaMjpsi.Rebin(5)
        print(eephoton12_data_deltaMjpsi.GetNbinsX())
        eephoton12_data_deltaMjpsi.Scale(1, "width")
    
    if option == 33:
        data_dileptonphoton = list_data2.FindObject("DileptonPhotonInvMass_cut")
        dileptonphoton_deltaM_jpsi = data_dileptonphoton.FindObject("DeltaMass_Jpsi");
        dileptonphoton_deltaM_jpsi.SetDirectory(0);
        make_common_style(dileptonphoton_deltaM_jpsi, kFullCrossX, 1.0, kPink+7, 1, 0);
        ROOT.SetOwnership(dileptonphoton_deltaM_jpsi, False);
        dileptonphoton_deltaM_jpsi.Rebin(5)
        print(dileptonphoton_deltaM_jpsi.GetNbinsX())
        dileptonphoton_deltaM_jpsi.Scale(1, "width")
        data_eephoton_chic12 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic12")
        eephoton12_data_deltaMjpsi = data_eephoton_chic12.FindObject("DeltaMass_Jpsi");
        eephoton12_data_deltaMjpsi.SetDirectory(0);
        make_common_style(eephoton12_data_deltaMjpsi, kFullCrossX, 1.0, kMagenta+3, 1, 0);
        ROOT.SetOwnership(eephoton12_data_deltaMjpsi, False);
        eephoton12_data_deltaMjpsi.Rebin(5)
        print(eephoton12_data_deltaMjpsi.GetNbinsX())
        eephoton12_data_deltaMjpsi.Scale(1, "width")
    
    if option == 50: 
        data_eephoton_chic12 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic12") 
        eephoton12_data_M = data_eephoton_chic12.FindObject("Mass_DiletonPhoton");
        eephoton12_data_M.SetDirectory(0);
        make_common_style(eephoton12_data_M, kFullCrossX, 1.0, kBlue+2, 1, 0);
        ROOT.SetOwnership(eephoton12_data_M, False);
        eephoton12_data_M.Rebin(5)
        print(eephoton12_data_M.GetNbinsX())
        eephoton12_data_M.Scale(1, "width")
        data_eephoton_chic12 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic12")
        eephoton12_data_deltaMjpsi = data_eephoton_chic12.FindObject("DeltaMass_Jpsi");
        eephoton12_data_deltaMjpsi.SetDirectory(0);
        make_common_style(eephoton12_data_deltaMjpsi, kFullCrossX, 1.0, kMagenta+3, 1, 0);
        ROOT.SetOwnership(eephoton12_data_deltaMjpsi, False);
        eephoton12_data_deltaMjpsi.Rebin(5)
        print(eephoton12_data_deltaMjpsi.GetNbinsX())
        eephoton12_data_deltaMjpsi.Scale(1, "width")


   
    ymin = 1
    ymax = 100
    if option == 0:
        ymax =eephoton12_deltaM.GetMaximum()
    if option == 1:
        ymax =eephoton2_deltaM.GetMaximum()
    if (option == 2) or (option == 4) or (option == 5): 
        ymax =eephoton12_deltaM.GetMaximum()
    if (option == 3 or option == 30 ):
        ymax = dileptonphoton_deltaM.GetMaximum()  
    if option == 12:
        ymax = eephoton12_data_deltaM.GetMaximum()
    if option == 20:
        ymax =eephoton12_deltaM.GetMaximum()
    if option == 23: 
        ymax = eephoton12_deltaMjpsi.GetMaximum()*2
        ymin = 250000
    if option == 26:
        ymax = eephoton12_data_deltaMjpsi.GetMaximum()
        ymin = 0
    if ( option == 33):
        ymax = dileptonphoton_deltaM_jpsi.GetMaximum()*5
        ymin = 999
    if option ==13:
        ymin = 89999
        ymax = dileptonphoton_deltaM_jpsi.GetMaximum()*1.5
    if option == 50: 
        ymax = max(eephoton12_data_deltaMjpsi.GetMaximum(), eephoton12_data_M.GetMaximum())
        ymin = 0
    c1 = TCanvas("DeltaMass","\Delta Mass",0,0,900,900);
    p1 = c1.cd();
    p1.SetPad(0,0.01,0.99,1);
    p1.SetMargin(0.15,0.05,0.12,0.05);
    p1.SetTicks(1,1);
    if (option != 50 and option != 26): 
        p1.SetLogy() #log scale
    xmin = 0
    xmax = 1
    if option ==12:
        xmin = 0.3 
        xmax = 0.6
    if (option == 23  ):
        xmin = 3.46
        xmax = 3.61
    if option == 13: 
        xmin = 3.2
        xmax = 4
    if ( option == 33 ): 
        xmin = 3.4
        xmax = 3.6
    if (option == 26 ): 
        xmin = 3.4
        xmax = 3.6
    if option == 50:
        xmin =3.4
        xmax = 3.6

    frame1 = p1.DrawFrame(xmin, ymin ,xmax, ymax+(ymax*0.2));
    if (option == 23 or option == 26 or option == 13 or option == 33): 
        frame1.GetXaxis().SetTitle("\Deltam + m_{J/\psi}^{PDG} [GeV/c^{2}]");
    elif option == 50:
        frame1.GetXaxis().SetTitle("m_{\gamma e^{+} e^{-}} [GeV/c^{2}]");
    else:
        frame1.GetXaxis().SetTitle("\DeltaM [GeV/c^{2}]");
    frame1.GetYaxis().SetTitle("Counts per GeV/c^{2}");
    frame1.GetXaxis().SetTitleSize(0.045);
    frame1.GetYaxis().SetTitleSize(0.045);
    frame1.GetXaxis().SetTitleOffset(1.1);
    frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetMaxDigits(3);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    ROOT.SetOwnership(frame1,False);


    if option == 0: 
        leg = TLegend(0.6,0.6,1.0,0.75);
        eephoton1_deltaM.Draw("Esame")
        leg.AddEntry(eephoton1_deltaM, "\gamma e^{+} e^{-} from \chi_{c1} MC truth","LP");
    if option == 1: 
        leg = TLegend(0.6,0.6,1.0,0.75);
        eephoton2_deltaM.Draw("Esame")
        leg.AddEntry(eephoton2_deltaM, "\gamma e^{+} e^{-} from \chi_{c2} MC truth","LP");
    if option == 2: 
        leg = TLegend(0.6,0.6,1.,0.7);
        eephoton12_deltaM.Draw("Esame")
        eephoton12_deltaM.Draw("Esame,hist")
        eephoton1_deltaM.Draw("Esame")
        eephoton2_deltaM.Draw("Esame")
        leg.AddEntry(eephoton12_deltaM, "\gamma e^{+} e^{-} from \chi_{c1} and \chi_{c2} ","LP")
        leg.AddEntry(eephoton1_deltaM, "\gamma e^{+} e^{-} from \chi_{c1} ","LP");
        leg.AddEntry(eephoton2_deltaM, "\gamma e^{+} e^{-} from \chi_{c2} ","LP");
    if option == 4: 
        leg = TLegend(0.6,0.65,1.0,0.75);
        eephoton1_deltaM.Draw("Esame")
        eephoton2_deltaM.Draw("Esame")
        eephoton12_deltaM.Draw("Esame")
        leg.AddEntry(eephoton1_deltaM, "\gamma e^{+} e^{-} from \chi_{c1} MC truth","LP");
        leg.AddEntry(eephoton2_deltaM, "\gamma e^{+} e^{-} from \chi_{c2} MC truth","LP");
        leg.AddEntry(eephoton12_deltaM, "\gamma e^{+} e^{-} from \chi_{c} MC truth","LP")
    if option == 3: 
        leg = TLegend(0.6,0.6,1.0,0.75);
        dileptonphoton_deltaM.Draw("Esame")
        leg.AddEntry(dileptonphoton_deltaM, "\gamma e^{+} e^{-}","LP")
    if option == 5: 
        leg = TLegend(0.6,0.6,1.0,0.75);
        eephoton1_deltaM.Draw("Esame")
        leg.AddEntry(eephoton1_deltaM, "\gamma e^{+} e^{-} from \chi_{c1}","LP");
        eephoton2_deltaM.Draw("Esame")
        leg.AddEntry(eephoton2_deltaM, "\gamma e^{+} e^{-} from \chi_{c2}","LP");

    if option == 12: 
        leg = TLegend(0.6,0.6,0.95,0.7);
        eephoton12_data_deltaM.Draw("Esame")
        eephoton12_data_deltaM.Draw("Esame,hist")
        eephoton1_data_deltaM.Draw("Esame")
        eephoton2_data_deltaM.Draw("Esame")
        leg.AddEntry(eephoton12_data_deltaM, "\gamma e^{+} e^{-} from \chi_{c1} and \chi_{c2} ","LP")
        leg.AddEntry(eephoton1_data_deltaM, "\gamma e^{+} e^{-} from \chi_{c1} ","LP");
        leg.AddEntry(eephoton2_data_deltaM, "\gamma e^{+} e^{-} from \chi_{c2} ","LP");

    if option == 13: 
        leg = TLegend(0.7,0.45,1.0,0.5);
        dileptonphoton_deltaM_jpsi.Draw("Esame")
        leg.AddEntry(dileptonphoton_deltaM_jpsi, "\gamma e^{+} e^{-}","LP")
    
    if option == 21: 
        leg = TLegend(0.6,0.6,1.0,0.75);
        eephoton1_deltaMjpsi.Draw("Esame")
        leg.AddEntry(eephoton1_deltaMjpsi, "\gamma e^{+} e^{-} from \chi_{c1}","LP");
        eephoton2_deltaMjpsi.Draw("Esame")
        leg.AddEntry(eephoton2_deltaMjpsi, "\gamma e^{+} e^{-} from \chi_{c2}","LP");
    if option == 22: 
        leg = TLegend(0.35,0.2,0.75,0.3);
        eephoton12_deltaMjpsi.Draw("Esame")
        eephoton12_deltaMjpsi.Draw("Esame,hist")
        eephoton1_deltaMjpsi.Draw("Esame")
        eephoton2_deltaMjpsi.Draw("Esame")
        leg.AddEntry(eephoton12_deltaMjpsi, "\gamma e^{+} e^{-} from \chi_{c1} and \chi_{c2} ","LP")
        leg.AddEntry(eephoton1_deltaMjpsi, "\gamma e^{+} e^{-} from \chi_{c1} ","LP");
        leg.AddEntry(eephoton2_deltaMjpsi, "\gamma e^{+} e^{-} from \chi_{c2} ","LP");
    if option == 23: 
        leg = TLegend(0.17,0.825,0.5,0.925);
        eephoton12_deltaMjpsi.Draw("Esame")
        eephoton12_deltaMjpsi.Draw("Esame,hist")
        eephoton1_deltaMjpsi.Draw("Esame")
        eephoton2_deltaMjpsi.Draw("Esame")
        leg.AddEntry(eephoton12_deltaMjpsi, "\chi_{c1,c2} \\rightarrow \gamma e^{+} e^{-}","LP")
        leg.AddEntry(eephoton1_deltaMjpsi, "\chi_{c1} \\rightarrow \gamma e^{+} e^{-}","LP");
        leg.AddEntry(eephoton2_deltaMjpsi, "\chi_{c2} \\rightarrow \gamma e^{+} e^{-}","LP");

    if option == 26: 
        leg = TLegend(0.17,0.625,0.5,0.725);
        eephoton12_data_deltaMjpsi.Draw("Esame")
        eephoton12_data_deltaMjpsi.Draw("Esame,hist")
        eephoton1_data_deltaMjpsi.Draw("Esame")
        eephoton2_data_deltaMjpsi.Draw("Esame")
        leg.AddEntry(eephoton12_data_deltaMjpsi, "\chi_{c1,c2} \\rightarrow \gamma e^{+} e^{-}","LP")
        leg.AddEntry(eephoton1_data_deltaMjpsi, "\chi_{c1} \\rightarrow \gamma e^{+} e^{-}","LP");
        leg.AddEntry(eephoton2_data_deltaMjpsi, "\chi_{c2} \\rightarrow \gamma e^{+} e^{-}","LP");
    
    if option == 30: 
        leg = TLegend(0.6,0.6,1.0,0.75);
        dileptonphoton_deltaM.Draw("Esame")
        eephoton12_data_deltaM.Draw("Esame")
        leg.AddEntry(dileptonphoton_deltaM, "\gamma e^{+} e^{-}","LP")
        leg.AddEntry(eephoton12_data_deltaM, "\gamma e^{+} e^{-} from \chi_{c1} and \chi_{c2} matched MC","LP")

    if option == 33: 
        leg = TLegend(0.18,0.825,0.5,0.875);
        dileptonphoton_deltaM_jpsi.Draw("Esame")
        eephoton12_data_deltaMjpsi.Draw("Esame,hist")
        leg.AddEntry(dileptonphoton_deltaM_jpsi, "\gamma e^{+} e^{-}","LP")
        leg.AddEntry(eephoton12_data_deltaMjpsi, "\chi_{c1,c2} \\rightarrow \gamma e^{+} e^{-} matched MC","LP")
    
    if option == 50: 
        leg = TLegend(0.17, 0.8,0.5,0.9);
        eephoton12_data_deltaMjpsi.Draw("Esame,hist")
        eephoton12_data_M.Draw("Esame, hist")
        leg.AddEntry(eephoton12_data_deltaMjpsi, "\Deltam_{\gamma e^{+} e^{-}}^{\chi_{c}} + m_{J/\psi}^{PDG} ","LP")
        leg.AddEntry(eephoton12_data_M, "m_{\gamma e^{+} e^{-}}^{\chi_{c}}","LP")
    

    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    if ( option == 33 or option == 26 ):
        leg.SetTextSize(0.025);
    elif option == 23: 
        leg.SetTextSize(0.025);
    else:
        leg.SetTextSize(0.03);
    leg.Draw("");
    ROOT.SetOwnership(leg,False);   
    if (option == 13):
        txt = TPaveText(0.90,0.85,0.9,0.95,"NDC"); 
    # elif (option == 33):
    #     txt = TPaveText(0.90,0.7,0.9,0.9,"NDC"); 
    elif ( option == 26):
        txt = TPaveText(0.4,0.85,0.4,0.95,"NDC");
    else:
        txt = TPaveText(0.90,0.85,0.9,0.95,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(33);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.03);
    txt.AddText("ALICE simulation");
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    if (option == 13):
        txt2 = TPaveText(0.845,0.8,0.83,0.925,"NDC");
    # elif (option == 33):
    #     txt2 = TPaveText(0.845,0.7,0.83,0.825,"NDC");
    elif (option == 26):
        txt2 = TPaveText(0.345,0.8,0.33,0.925,"NDC");
    else: 
        txt2 = TPaveText(0.845,0.8,0.83,0.925,"NDC");
    txt2.SetFillColor(kWhite);
    txt2.SetFillStyle(0);
    txt2.SetBorderSize(0);
    txt2.SetTextAlign(33);#middle,left
    txt2.SetTextFont(42);#helvetica
    txt2.SetTextSize(0.03);
    txt2.AddText("this thesis");
    txt2.Draw();
    ROOT.SetOwnership(txt2,False);
    if ( option == 13):
        txt3 = TPaveText(0.9,0.75,0.9,0.90,"NDC");
    # elif ( option == 33):
    #     txt3 = TPaveText(0.9,0.65,0.9,0.8,"NDC");
    elif (  option == 26):
        txt3 = TPaveText(0.4,0.75,0.4,0.9,"NDC");
    else: 
        txt3 = TPaveText(0.9,0.75,0.9,0.90,"NDC");
    txt3.SetFillColor(kWhite);
    txt3.SetFillStyle(0);
    txt3.SetBorderSize(0);
    txt3.SetTextAlign(33);#middle,left
    txt3.SetTextFont(42);#helvetica
    txt3.SetTextSize(0.03);
    txt3.AddText("pp, #sqrt{s} = 13.6TeV");
    txt3.Draw();
    ROOT.SetOwnership(txt3,False);
    if (option == 12 or option == 50): 
        txt4 = TPaveText(0.8,0.77,0.87,0.8,"NDC");
        txt4.SetFillColor(kWhite);
        txt4.SetFillStyle(0);
        txt4.SetBorderSize(0);
        txt4.SetTextAlign(33);#middle,left
        txt4.SetTextFont(42);#helvetica
        txt4.SetTextSize(0.03);
        txt4.AddText("MC matched");
        txt4.Draw();
        ROOT.SetOwnership(txt4,False);
    if ( option == 26): 
        txt4 = TPaveText(0.3,0.77,0.37,0.8,"NDC");
        txt4.SetFillColor(kWhite);
        txt4.SetFillStyle(0);
        txt4.SetBorderSize(0);
        txt4.SetTextAlign(33);#middle,left
        txt4.SetTextFont(42);#helvetica
        txt4.SetTextSize(0.03);
        txt4.AddText("MC matched");
        txt4.Draw();
        ROOT.SetOwnership(txt4,False);
    if (option == 2 or option == 23): 
        txt4 = TPaveText(0.88,0.77,0.88,0.8,"NDC");
        txt4.SetFillColor(kWhite);
        txt4.SetFillStyle(0);
        txt4.SetBorderSize(0);
        txt4.SetTextAlign(33);#middle,left
        txt4.SetTextFont(42);#helvetica
        txt4.SetTextSize(0.03);
        txt4.AddText("MC generated");
        txt4.Draw();
        ROOT.SetOwnership(txt4,False);
    if option == 23:
        arrow = TArrow( 3.51069, 600000, 3.51069, 250000, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
        arrow2 = TArrow( 3.55617, 600000, 3.55617, 250000, 0.02, '|>' )
        arrow2.SetFillStyle( 1001 )
        arrow2.Draw()
    if option == 26:
        arrow = TArrow( 3.51069, 10000, 3.51069, 0, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
        arrow2 = TArrow( 3.55617, 10000, 3.55617, 0, 0.02, '|>' )
        arrow2.SetFillStyle( 1001 )
        arrow2.Draw()
    elif option == 13:
        arrow = TArrow( 3.51069, 120000, 3.51069, 89999, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
        arrow2 = TArrow( 3.55617, 120000, 3.55617, 89999, 0.02, '|>' )
        arrow2.SetFillStyle( 1001 )
        arrow2.Draw()
    elif option == 33:
        arrow = TArrow( 3.51069, 2000, 3.51069, 999, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
        arrow2 = TArrow( 3.55617, 2000, 3.55617, 999, 0.02, '|>' )
        arrow2.SetFillStyle( 1001 )
        arrow2.Draw()
    elif option == 50: 
        arrow = TArrow( 3.51069, 10000, 3.51069, 10, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
        arrow2 = TArrow( 3.55617, 10000, 3.55617, 10, 0.02, '|>' )
        arrow2.SetFillStyle( 1001 )
        arrow2.Draw()
    

    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs(plotname);





if __name__ == "__main__":
    filename = "AnalysisResults_chicall_20240224.root"
    deltamass_plot(filename, 23, "20240225/plot_deltamass_eegammachic12_truth.pdf")
    deltamass_plot(filename, 23, "20240225/plot_deltamass_eegammachic12_truth.svg")
    deltamass_plot(filename, 26, "20240225/plot_deltamass_eegammachic12_matched.pdf")
    deltamass_plot(filename, 26, "20240225/plot_deltamass_eegammachic12_matched.svg")
    deltamass_plot(filename, 13, "20240225/plot_deltamass_eegamma_triple.pdf")
    deltamass_plot(filename, 13, "20240225/plot_deltamass_eegamma_triple.svg")
    deltamass_plot(filename, 33, "20240225/plot_deltamass_eegammachic12_triplematched.pdf")
    deltamass_plot(filename, 33, "20240225/plot_deltamass_eegammachic12_triplematched.svg")
    deltamass_plot(filename, 50, "20240225/plot_deltamass_mass_eegammachic12_matched.pdf")
    deltamass_plot(filename, 50, "20240225/plot_deltamass_mass_eegammachic12_matched.svg")
    
# option 2: delta Mass truth MC for chic1, chic2 and chic12 together
# option 12: delta Mass matched MC for chic1, chic2 and chic12 together 