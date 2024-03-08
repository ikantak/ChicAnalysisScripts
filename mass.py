import re
import numpy as np
import datetime
import math

import ROOT
import ctypes
from ctypes import *
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency, TArrow
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kGray, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
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

def mass_plot(filename1, option, plotname):
    rootfile_data = TFile.Open(filename1, "READ");
    list_data = rootfile_data.Get("analysis-dilepton-photon");
    list_data2 = list_data.Get("output");
    list_same = rootfile_data.Get("analysis-same-event-pairing")
    list_same2 = list_same.Get("output")

    if (option == 0 or option == 11): 
        data_cut_ee_jpsi = list_data2.FindObject("DileptonsSelected_cut_matchedMC_eeFromJpsi")
        eejpsi_cut_data_mass = data_cut_ee_jpsi.FindObject("Mass_dilepton")
        eejpsi_cut_data_mass.SetDirectory(0)
        make_common_style(eejpsi_cut_data_mass, kFullCross, 1.0, kCyan+2, 1, 0)
        ROOT.SetOwnership(eejpsi_cut_data_mass, False)
        eejpsi_cut_data_mass.Rebin(5)
        print(eejpsi_cut_data_mass.GetNbinsX())
        eejpsi_cut_data_mass.Scale(1, "width")

    if (option == 1 or option == 5): 
       
        data_cut_eephoton_chic1 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic1")
        data_cut_eephoton_chic2 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic2")

        eephoton1_cut_data_mass = data_cut_eephoton_chic1.FindObject("Mass_DiletonPhoton");
        eephoton1_cut_data_mass.SetDirectory(0);
        make_common_style(eephoton1_cut_data_mass, kFullCrossX, 1.0, kViolet-1, 1, 0);
        ROOT.SetOwnership(eephoton1_cut_data_mass, False);
        print(eephoton1_cut_data_mass.GetNbinsX())
        eephoton1_cut_data_mass.Rebin(10)
        print(eephoton1_cut_data_mass.GetNbinsX())
        eephoton1_cut_data_mass.Scale(1, "width")

        eephoton2_cut_data_mass = data_cut_eephoton_chic2.FindObject("Mass_DiletonPhoton");
        eephoton2_cut_data_mass.SetDirectory(0);
        make_common_style(eephoton2_cut_data_mass, kFullCrossX, 1.0, kViolet-4, 1, 0);
        ROOT.SetOwnership(eephoton2_cut_data_mass, False)
        print(eephoton2_cut_data_mass.GetNbinsX())
        eephoton2_cut_data_mass.Rebin(10)
        print(eephoton2_cut_data_mass.GetNbinsX())
        eephoton2_cut_data_mass.Scale(1, "width")

    
    if (option == 2 or option == 4 or option == 5): 
        data_cut_eephoton_chic12 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic12")
        eephoton12_cut_data_mass = data_cut_eephoton_chic12.FindObject("Mass_DiletonPhoton");
        make_common_style(eephoton12_cut_data_mass, kFullCrossX, 1.0, kMagenta+3, 1, 0);
        ROOT.SetOwnership(eephoton12_cut_data_mass, False);
        print(eephoton12_cut_data_mass.GetNbinsX())
        eephoton12_cut_data_mass.Rebin(10)
        print(eephoton12_cut_data_mass.GetNbinsX())
        #binwidth = eephoton12_cut_data_mass.GetXaxis().GetBinWidth(1)
        eephoton12_cut_data_mass.Scale(1, "width")

    if (option == 3 or option == 4): 
        data_cut_eephoton = list_data2.FindObject("DileptonPhotonInvMass_cut")
        eephoton_cut_data_mass = data_cut_eephoton.FindObject("Mass_DiletonPhoton");
        eephoton_cut_data_mass.SetDirectory(0);
        make_common_style(eephoton_cut_data_mass, kFullCrossX, 1.0, kPink+7, 1, 0);
        ROOT.SetOwnership(eephoton_cut_data_mass, False);
        print(eephoton_cut_data_mass.GetNbinsX())
        eephoton_cut_data_mass.Rebin(10)
        print(eephoton_cut_data_mass.GetNbinsX())
        eephoton_cut_data_mass.Scale(1, "width")
        
        #eephoton_cut_data_mass.Rebin(10)
    if (option == 10 or option == 11): 
        data_cut_ee_jpsi = list_data2.FindObject("MCTruthGenPair_eeFromJpsi")
        eejpsi_cut_mass = data_cut_ee_jpsi.FindObject("Mass")
        eejpsi_cut_mass.SetDirectory(0)
        make_common_style(eejpsi_cut_mass, kFullCross, 1.0, kBlue+2, 1, 0)
        ROOT.SetOwnership(eejpsi_cut_mass, False)
        eejpsi_cut_mass.Rebin(5)
        print(eejpsi_cut_mass.GetNbinsX())
        eejpsi_cut_mass.Scale(1, "width")

    if option == 12: 
        data_cut_ee_chic1 = list_data2.FindObject("DileptonsSelected_cut_matchedMC_eeFromChic1")
        data_cut_ee_chic2 = list_data2.FindObject("DileptonsSelected_cut_matchedMC_eeFromChic2")

        ee1_cut_data_mass = data_cut_ee_chic1.FindObject("Mass");
        ee1_cut_data_mass.SetDirectory(0);
        make_common_style(ee1_cut_data_mass, kFullCross, 1.0, kViolet-1, 1, 0);
        ROOT.SetOwnership(ee1_cut_data_mass, False);
        #print(eephoton1_cut_data_mass.GetNbinsX())
        ee1_cut_data_mass.Scale(1, "width")

        ee2_cut_data_mass = data_cut_ee_chic2.FindObject("Mass");
        ee2_cut_data_mass.SetDirectory(0);
        make_common_style(ee2_cut_data_mass, kFullCross, 1.0, kViolet-4, 1, 0);
        ROOT.SetOwnership(ee2_cut_data_mass, False)
        ee2_cut_data_mass.Scale(1, "width")

    if option == 13:
        mc_cut_ee_chic1 = list_data2.FindObject("MCTruthGenPair_cut_matchedMC_eeFromChic1")
        mc_cut_ee_chic2 = list_data2.FindObject("MCTruthGenPair_cut_matchedMC_eeFromChic2")

        ee1_cut_mc_mass = mc_cut_ee_chic1.FindObject("Mass");
        ee1_cut_mc_mass.SetDirectory(0);
        make_common_style(ee1_cut_mc_mass, kFullCross, 1.0, kViolet-1, 1, 0);
        ROOT.SetOwnership(ee1_cut_mc_mass, False);
        #print(eephoton1_cut_data_mass.GetNbinsX())
        ee1_cut_mc_mass.Scale(1, "width")

        ee2_cut_mc_mass = mc_cut_ee_chic2.FindObject("Mass");
        ee2_cut_mc_mass.SetDirectory(0);
        make_common_style(ee2_cut_mc_mass, kFullCross, 1.0, kViolet-4, 1, 0);
        ROOT.SetOwnership(ee2_cut_mc_mass, False)
        ee2_cut_mc_mass.Scale(1, "width") 

    if option == 15: 
        mc_cut_ee = list_data2.FindObject("DileptonsSelected_cut")
        mc_cut_ee_jpsi = list_data2.FindObject("DileptonsSelected_cut_matchedMC_eeFromJpsi")

        ee_cut_mass = mc_cut_ee.FindObject("Mass_dilepton");
        ee_cut_mass.SetDirectory(0);
        make_common_style(ee_cut_mass, kFullCross, 1.0, kPink+8, 1, 0);
        ROOT.SetOwnership(ee_cut_mass, False);
        ee_cut_mass.Rebin(5)
        print(ee_cut_mass.GetNbinsX())
        #print(eephoton1_cut_data_mass.GetNbinsX())
        ee_cut_mass.Scale(1, "width")

        eejpsi_cut_mass = mc_cut_ee_jpsi.FindObject("Mass_dilepton");
        eejpsi_cut_mass.SetDirectory(0);
        make_common_style(eejpsi_cut_mass, kFullCross, 1.0, kCyan+3, 1, 0);
        ROOT.SetOwnership(eejpsi_cut_mass, False)
        eejpsi_cut_mass.Rebin(5)
        print(eejpsi_cut_mass.GetNbinsX())
        eejpsi_cut_mass.Scale(1, "width")

    if option == 20 or option == 22: 
        pm_ee = list_same2.FindObject("PairsBarrelSEPM_jpsiO2MCdebugCuts2")
        pp_ee = list_same2.FindObject("PairsBarrelSEPP_jpsiO2MCdebugCuts2")
        mm_ee = list_same2.FindObject("PairsBarrelSEMM_jpsiO2MCdebugCuts2")
        mc_cut_ee = list_data2.FindObject("DileptonsSelected_cut")

        ee_cut_mass = mc_cut_ee.FindObject("Mass_dilepton");
        ee_cut_mass.SetDirectory(0);
        make_common_style(ee_cut_mass, kFullCross, 1.0, kPink+6, 1, 0);
        ROOT.SetOwnership(ee_cut_mass, False);
        ee_cut_mass.Rebin(5)
        print(ee_cut_mass.GetNbinsX())
        #print(eephoton1_cut_data_mass.GetNbinsX())
        ee_cut_mass.Scale(1, "width")

        ee_pm_mass = pm_ee.FindObject("Mass_dilepton");
        ee_pm_mass.SetDirectory(0);
        make_common_style(ee_pm_mass, kFullCross, 1.0, kPink+7, 1, 0);
        ROOT.SetOwnership(ee_pm_mass, False);
        ee_pm_mass.Rebin(5)
        print(ee_pm_mass.GetNbinsX())
        ee_pm_mass.Scale(1, "width")
        
        ee_pp_mass = pp_ee.FindObject("Mass_dilepton");
        ee_mm_mass = mm_ee.FindObject("Mass_dilepton");
        background = ee_pp_mass.Clone("background")
        background.Add(ee_mm_mass)
        make_common_style(background, kFullCross, 1.0, kBlack, 1, 0);
        ROOT.SetOwnership(background, False)
        background.Rebin(5)
        background.Scale(1, "width")
    
    if option == 21 or option == 25: 
        pm_ee = list_same2.FindObject("PairsBarrelSEPM_jpsiO2MCdebugCuts2")
        pp_ee = list_same2.FindObject("PairsBarrelSEPP_jpsiO2MCdebugCuts2")
        mm_ee = list_same2.FindObject("PairsBarrelSEMM_jpsiO2MCdebugCuts2")

        ee_pm_mass = pm_ee.FindObject("Mass_dilepton");
        ee_pm_mass.SetDirectory(0);
        make_common_style(ee_pm_mass, kFullCross, 1.0, kPink+8, 1, 0);
        ROOT.SetOwnership(ee_pm_mass, False);
        #print(eephoton1_cut_data_mass.GetNbinsX())
        
        ee_pp_mass = pp_ee.FindObject("Mass_dilepton");
        ee_mm_mass = mm_ee.FindObject("Mass_dilepton");
        background = ee_pp_mass.Clone("background")
        background.Add(ee_mm_mass, 1)
        make_common_style(background, kFullCross, 1.0, kCyan+3, 1, 0);
        ROOT.SetOwnership(background, False)
        ee_withoutbackground = ee_pm_mass.Clone("ee without background")
        ee_withoutbackground.Add(background, -1)
        ee_withoutbackground.Rebin(5)
        print(ee_withoutbackground.GetNbinsX())
        ee_withoutbackground.Scale(1, "width")

        mc_cut_ee_jpsi = list_data2.FindObject("DileptonsSelected_cut_matchedMC_eeFromJpsi")
        eejpsi_cut_mass = mc_cut_ee_jpsi.FindObject("Mass_dilepton");
        eejpsi_cut_mass.SetDirectory(0);
        make_common_style(eejpsi_cut_mass, kFullCross, 1.0, kCyan+3, 1, 0);
        ROOT.SetOwnership(eejpsi_cut_mass, False)
        eejpsi_cut_mass.Rebin(5)
        print(eejpsi_cut_mass.GetNbinsX())
        eejpsi_cut_mass.Scale(1, "width")
    
    if option == 30: 
        mc_cut_eephoton_chic1 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic1")
        mc_cut_eephoton_chic2 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic2")

        eephoton1_cut_mc_mass = mc_cut_eephoton_chic1.FindObject("Mass_DileptonPhoton");
        eephoton1_cut_mc_mass.SetDirectory(0);
        make_common_style(eephoton1_cut_mc_mass, kFullCrossX, 1.0, kRed, 1, 0);
        ROOT.SetOwnership(eephoton1_cut_mc_mass, False);
        print(eephoton1_cut_mc_mass.GetNbinsX())
        #eephoton1_cut_mc_mass.Rebin(5)
        print(eephoton1_cut_mc_mass.GetNbinsX())
        eephoton1_cut_mc_mass.Scale(1, "width")

        eephoton2_cut_mc_mass = mc_cut_eephoton_chic2.FindObject("Mass_DileptonPhoton");
        eephoton2_cut_mc_mass.SetDirectory(0);
        make_common_style(eephoton2_cut_mc_mass, kFullCrossX, 1.0, kOrange+7, 1, 0);
        ROOT.SetOwnership(eephoton2_cut_mc_mass, False)
        print(eephoton2_cut_mc_mass.GetNbinsX())
        #eephoton2_cut_mc_mass.Rebin(5)
        print(eephoton2_cut_mc_mass.GetNbinsX())
        eephoton2_cut_mc_mass.Scale(1, "width")

        mc_cut_eephoton_chic12 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic12")
        eephoton12_cut_mc_mass = mc_cut_eephoton_chic12.FindObject("Mass_DileptonPhoton");
        make_common_style(eephoton12_cut_mc_mass, kFullCrossX, 1.0, kRed+2, 1, 0);
        ROOT.SetOwnership(eephoton12_cut_mc_mass, False);
        print(eephoton12_cut_mc_mass.GetNbinsX())
        #eephoton12_cut_mc_mass.Rebin(5)
        print(eephoton12_cut_mc_mass.GetNbinsX())
        eephoton12_cut_mc_mass.Scale(1, "width")

    if option == 11: 
        #h1eff = TH1D("mass", "mass", 4500, 0, 4.50)
        # eejpsi_cut_data_mass.Rebin(5)
        # eejpsi_cut_mass.Rebin(5)
        max1= eejpsi_cut_mass.GetMaximum()
        max2 = eejpsi_cut_data_mass.GetMaximum()
        eejpsi_cut_mass.Scale(1./max1)
        eejpsi_cut_data_mass.Scale(1./max2)
        # outfile = TFile("eejpsi_truthmatched.root", "RECREATE");
        # outfile.WriteTObject(eejpsi_cut_mass);
        # outfile.WriteTObject(eejpsi_cut_data_mass)
        # outfile.Close();
    if option == 25:
        ratio = eejpsi_cut_mass.Clone()
        ratio.Sumw2()
        ratio.Divide(ee_withoutbackground, eejpsi_cut_mass, 1, 1, option = "B")
        make_common_style(ratio, kFullCrossX, 1.0, kBlack, 1, 0);
        ROOT.SetOwnership(ratio, False);

    if option == 0: 
        ymax = eejpsi_cut_data_mass.GetMaximum()*1.1
        ymin = eejpsi_cut_data_mass.GetMinimum()
    if option == 1: 
        ymax = max(eephoton1_cut_data_mass.GetMaximum(), eephoton2_cut_data_mass.GetMaximum())
        ymin = min(eephoton1_cut_data_mass.GetMinimum(), eephoton2_cut_data_mass.GetMinimum())
    if option == 2: 
        ymax = eephoton12_cut_data_mass.GetMaximum()
        ymin = eephoton12_cut_data_mass.GetMinimum()
    if option == 3:
        ymax = eephoton_cut_data_mass.GetMaximum()*1.25
        ymin = 300000
    if option == 4: 
        ymax = eephoton_cut_data_mass.GetMaximum()*5
        ymin = 800
    if option == 5:
        ymax = eephoton12_cut_data_mass.GetMaximum()+3500
        ymin = min(eephoton1_cut_data_mass.GetMinimum(), eephoton2_cut_data_mass.GetMinimum())
    if option == 10: 
        ymax = eejpsi_cut_mass.GetMaximum()*2
        ymin = 99999
    if option == 11: 
        ymax =1
        ymin = 0.0005
    if option == 12: 
        ymax = max(ee1_cut_data_mass.GetMaximum(), ee2_cut_data_mass.GetMaximum())
        ymin =min(ee1_cut_data_mass.GetMinimum(), ee2_cut_data_mass.GetMinimum())

    if option ==13:
        ymax =max(ee1_cut_mc_mass.GetMaximum(), ee2_cut_mc_mass.GetMaximum())
        ymin =min(ee1_cut_mc_mass.GetMinimum(), ee2_cut_mc_mass.GetMinimum())
    if option ==15:
        ymax =max(ee_cut_mass.GetMaximum(), eejpsi_cut_mass.GetMaximum())
        ymin = min(ee_cut_mass.GetMinimum(), eejpsi_cut_mass.GetMinimum())
    if option== 20 :
        ymax = max(ee_pm_mass.GetMaximum(), background.GetMaximum())+1000000
        ymin = min(ee_pm_mass.GetMinimum(), background.GetMinimum())
    if option == 22:
        ymax = max(ee_pm_mass.GetMaximum(), background.GetMaximum())
        ymin = min(ee_pm_mass.GetMinimum(), background.GetMinimum()) +754000
    if option == 21:
        ymax = max(ee_withoutbackground.GetMaximum(), eejpsi_cut_mass.GetMaximum())+1000000
        ymin = min(ee_withoutbackground.GetMinimum(), eejpsi_cut_mass.GetMinimum()) + 50000
    if option == 25: 
        ymax = 1.5
        ymin = 0.7
    if option == 30: 
        ymax = eephoton12_cut_mc_mass.GetMaximum()*2
        ymin = 999999

    if option == 20 or option == 0 or option == 10 or option == 21: 
        c1 = TCanvas("pT_distribution","pT distribution",0,0,1500,900);
    else: 
        c1 = TCanvas("pT_distribution","pT distribution",0,0,900,900);
    p1 = c1.cd();
    p1.SetPad(0,0.01,1,1);
    if option == 20 or option == 0 or option == 10 or option == 21: 
        p1.SetMargin(0.1,0.05,0.12,0.05);
    else: 
        p1.SetMargin(0.15,0.05,0.12,0.05);
    p1.SetTicks(1,1);
    if option == 10 or option == 11 or option == 4 or option == 30 or option == 3: 
        p1.SetLogy()
    #p1.SetLogy() #log scale
    if ymax <10 and option != 25:
        ymax += 1
    if (ymax >9 and ymax <100):
        ymax += 10
    if (ymax >=100 and ymax < 1000):
        ymax += 30
    if (ymax >= 1000 and ymax < 10000):
        ymax += 200
    if (ymax >= 10000 and ymax < 100000 and option != 0):
        ymax += 1000
    if (ymax >= 100000 and ymax < 1000000):
        ymax += 200000
    if (ymax >= 1000000 and ymax < 10000000):
        ymax += 200000

    # if (ymin < 1 or (ymin >= 1 and ymin < 10)):
    #     ymin = 1
    # elif (ymin >= 10 and ymin < 100):
    #     ymin = 100
    # elif (ymin >= 100 and ymin < 1000):
    #     ymin = 1000
    # elif (ymin>= 1000 and ymin < 10000):
    #     ymin = 10000
    # elif (ymin >= 10000 and ymin < 100000):
    #     ymin = 10000
    # elif (ymin >= 100000 and ymin < 1000000):
    #     ymin = 100000
    # # else: 
    #     ymin = 1000000
    
    xmax = 4
    xmin = 2
    if (option == 0):
        xmin = 2.5
        xmax = 3.2
    if (option == 2):
        xmin = 3.2
        xmax = 3.8
    if (option == 1): 
        xmin = 3.4
        xmax = 3.6
    if option == 3: 
        xmin = 3.2
        xmax = 4
    if option == 4: 
        xmin = 3.4
        xmax = 3.6
    if option == 5: 
        xmin = 3.4
        xmax = 3.6
    if option == 10: 
        xmin = 2
        xmax = 3.15
    if option == 11 or option == 15 or option == 21 or option == 25:
        xmin = 2.5        
        xmax = 3.2
    if option == 12: 
        xmin = 2
        xmax = 3.15
    if option == 13:
        xmin = 2
        xmax = 3.15
    if option == 20 : 
        xmin = 0
        xmax = 4.5
    if  option == 22:
        xmin = 2.35
        xmax = 3.25
    if option == 30:
        xmin = 3.47
        xmax = 3.57

    
    
    frame1 = p1.DrawFrame(xmin, ymin, xmax, ymax);
    if option == 25: 
        frame1.GetYaxis().SetTitle("ratio \\frac{MC matched J/\psi \\rightarrow e^{+}e^{-}}{e^{+}e^{-}-(e^{+}e^{+}+e^{-}e^{-})}");
    else: 
        frame1.GetYaxis().SetTitle("Counts per GeV/c^{2}");
    frame1.GetXaxis().SetTitleSize(0.045);
    frame1.GetYaxis().SetTitleSize(0.045);
    frame1.GetXaxis().SetTitleOffset(1.1);
    if option == 20 or option == 0 or option == 10 or option == 21: 
        frame1.GetYaxis().SetTitleOffset(1);
    else: 
        frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetMaxDigits(3);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    # if option == 5: 
    #     eephoton12_cut_data_mass.SetNdivisions(510,"x")
    

    #h1data.SetMarkerStyle(kFullCircle)
    if option == 0:
        leg = TLegend(0.125,0.7,0.34,0.75);
        eejpsi_cut_data_mass.Draw("Esame,hist")
        leg.AddEntry(eejpsi_cut_data_mass, "J/\psi \\rightarrow e^{+} e^{-}", "LP")
        frame1.GetXaxis().SetTitle("m_{e^{+}e^{-}} [GeV/c^{2}]");
    if option == 1: 
        leg = TLegend(0.2,0.2,0.5,0.3);
        eephoton1_cut_data_mass.Draw("Esame,hist")
        eephoton2_cut_data_mass.Draw("Esame,hist")
        leg.AddEntry(eephoton1_cut_data_mass, "\chi_{c1} \\rightarrow \gamma e^{+} e^{-}", "LP")
        leg.AddEntry(eephoton2_cut_data_mass, "\chi_{c2} \\rightarrow \gamma e^{+} e^{-}", "LP")
        frame1.GetXaxis().SetTitle("m_{\gamma e^{+}e^{-}} [GeV/c^{2}]");
    if option == 2: 
        leg = TLegend(0.2,0.2,0.5,0.3);
        eephoton12_cut_data_mass.Draw("Esame")
        leg.AddEntry(eephoton12_cut_data_mass, "\chi_{c1,c2} \\rightarrow \gamma e^{+} e^{-}", "LP")
        frame1.GetXaxis().SetTitle("m_{\gamma e^{+}e^{-}} [GeV/c^{2}]");
    if option == 3: 
        leg = TLegend(0.7,0.7,1.0,0.75);
        eephoton_cut_data_mass.Draw("Esame")
        leg.AddEntry(eephoton_cut_data_mass, "\gamma e^{+} e^{-} ", "LP")
        frame1.GetXaxis().SetTitle("m_{\gamma e^{+}e^{-}} [GeV/c^{2}]");
    if option == 4: 
        leg = TLegend(0.17,0.825,0.47,0.91);
        eephoton_cut_data_mass.Draw("Esame")
        eephoton12_cut_data_mass.Draw("Esame, hist")
        leg.AddEntry(eephoton_cut_data_mass, "\gamma e^{+} e^{-}", "LP")
        leg.AddEntry(eephoton12_cut_data_mass, "\chi_{c1,c2} \\rightarrow \gamma e^{+} e^{-} MC matched", "LP")  
        frame1.GetXaxis().SetTitle("m_{\gamma e^{+}e^{-}} [GeV/c^{2}]");
    if option == 5:
        leg = TLegend(0.17,0.8,0.5,0.91);
        eephoton12_cut_data_mass.Draw("Esame,hist")
        eephoton1_cut_data_mass.Draw("Esame,hist")
        eephoton2_cut_data_mass.Draw("Esame,hist")
        leg.AddEntry(eephoton12_cut_data_mass, "\chi_{c1,c2} \\rightarrow \gamma e^{+} e^{-}", "LP")
        leg.AddEntry(eephoton1_cut_data_mass, "\chi_{c1} \\rightarrow \gamma e^{+} e^{-}", "LP")
        leg.AddEntry(eephoton2_cut_data_mass, "\chi_{c2} \\rightarrow \gamma e^{+} e^{-}", "LP")
        frame1.GetXaxis().SetTitle("m_{\gamma e^{+}e^{-}} [GeV/c^{2}]");
    if option == 10:
        leg = TLegend(0.13,0.7,0.3,0.75);
        eejpsi_cut_mass.Draw("Esame,hist")
        leg.AddEntry(eejpsi_cut_mass, "J/\psi \\rightarrow e^{+} e^{-}", "LP")
        frame1.GetXaxis().SetTitle("m_{e^{+}e^{-}} [GeV/c^{2}]");
    if option == 11: 
        leg = TLegend(0.2,0.5,0.5,0.55);
        
        eejpsi_cut_data_mass.Draw("Esame,hist")
        eejpsi_cut_mass.Draw("Esame,hist")
        leg.AddEntry(eejpsi_cut_data_mass, "J/\psi \\rightarrow e^{+} e^{-} MC matched", "LP")
        leg.AddEntry(eejpsi_cut_mass, "e^{+} e^{-} from J/\psi MC generated", "LP")
        frame1.GetXaxis().SetTitle("m_{e^{+}e^{-}} [GeV/c^{2}]");
    if option == 12: 
        leg = TLegend(0.2,0.7,0.5,0.75);
        ee1_cut_data_mass.Draw("Esame")
        ee2_cut_data_mass.Draw("Esame")
        leg.AddEntry(ee1_cut_data_mass, "e^{+} e^{-} from \chi_{c1} MC matched", "LP")
        leg.AddEntry(ee2_cut_data_mass, "e^{+} e^{-} from \chi_{c2} MC matched", "LP")
        frame1.GetXaxis().SetTitle("m_{e^{+}e^{-}} [GeV/c^{2}]");
    if option == 13: 
        leg = TLegend(0.2,0.7,0.5,0.75);
        ee1_cut_mc_mass.Draw("Esame")
        ee2_cut_mc_mass.Draw("Esame")
        leg.AddEntry(ee1_cut_mc_mass, "e^{+} e^{-} from \chi_{c1} MC generated", "LP")
        leg.AddEntry(ee2_cut_mc_mass, "e^{+} e^{-} from \chi_{c2} MC generated", "LP")
        frame1.GetXaxis().SetTitle("m_{e^{+}e^{-}} [GeV/c^{2}]");
    if option == 15: 
        leg = TLegend(0.2,0.65,0.5,0.75);
        ee_cut_mass.Draw("Esame")
        eejpsi_cut_mass.Draw("Esame")
        leg.AddEntry(ee_cut_mass, "e^{+} e^{-}", "LP")
        leg.AddEntry(eejpsi_cut_mass, "e^{+} e^{-} from J/\psi MC matched", "LP")
        frame1.GetXaxis().SetTitle("m_{e^{+}e^{-}} [GeV/c^{2}]");
    if option == 20:
        leg = TLegend(0.15,0.7,0.3,0.8);
        # ee_cut_mass.Draw("Esame,hist")
        ee_pm_mass.Draw("Esame")
        background.Draw("Esame")
        # leg.AddEntry(ee_pm_mass, "e^{+} e^{-} with mass cut", "LP")
        leg.AddEntry(ee_pm_mass, "e^{+} e^{-}", "LP")
        leg.AddEntry(background, "e^{+} e^{+} and e^{-} e^{-}", "LP")
        frame1.GetXaxis().SetTitle("m_{e^{+}e^{-}} [GeV/c^{2}]");
        line = TLine( 2.5, ymax, 2.5, ymin)
        line.SetLineStyle( 2 )
        line.SetLineColor(13)
        line.Draw()
        line2 = TLine( 3.3, ymax, 3.3, ymin)
        line2.SetLineStyle(2)
        line2.SetLineColor(13)
        line2.Draw()
    if option== 22:
        leg = TLegend(0.2,0.65,0.5,0.75);
        # ee_cut_mass.Draw("Esame,hist")
        ee_pm_mass.Draw("Esame")
        background.Draw("Esame")
        # leg.AddEntry(ee_pm_mass, "e^{+} e^{-} with mass cut", "LP")
        leg.AddEntry(ee_pm_mass, "e^{+} e^{-}", "LP")
        leg.AddEntry(background, "e^{+} e^{+} and e^{-} e^{-}", "LP")
        frame1.GetXaxis().SetTitle("m_{e^{+}e^{-}} [GeV/c^{2}]");
        # line = TLine( 2.5, ymax, 2.5, ymin)
        # line.SetLineStyle( 1001 )
        # line.Draw()
        # line2 = TLine( 3.3, ymax, 3.3, ymin)
        # line2.SetLineStyle(1001 )
        #line2.Draw()
    if option == 21: 
        leg = TLegend(0.13,0.65,0.35,0.75);
        ee_withoutbackground.Draw("Esame,hist")
        eejpsi_cut_mass.Draw("Esame,hist")
        leg.AddEntry(ee_withoutbackground, "e^{+} e^{-} minus e^{+} e^{+} and e^{-} e^{-}", "LP")
        leg.AddEntry(eejpsi_cut_mass, "J/\psi \\rightarrow e^{+} e^{-} MC matched", "LP")
        frame1.GetXaxis().SetTitle("m_{e^{+}e^{-}} [GeV/c^{2}]");

    if option == 25: 
        leg = TLegend(0.2,0.65,0.5,0.75);
        ratio.Draw("Esame,hist")
        frame1.GetXaxis().SetTitle("m_{e^{+}e^{-}} [GeV/c^{2}]");
    
    if option == 30: 
        leg = TLegend(0.65,0.825,0.9,.925);
        eephoton12_cut_mc_mass.Draw("Esame,hist")
        eephoton1_cut_mc_mass.Draw("Esame")
        eephoton2_cut_mc_mass.Draw("Esame")
        leg.AddEntry(eephoton12_cut_mc_mass, "\chi_{c1,c2} \\rightarrow \gamma e^{+} e^{-}", "LP")
        leg.AddEntry(eephoton1_cut_mc_mass, "\chi_{c1} \\rightarrow \gamma e^{+} e^{-}", "LP")
        leg.AddEntry(eephoton2_cut_mc_mass, "\chi_{c1} \\rightarrow \gamma e^{+} e^{-}", "LP")
        frame1.GetXaxis().SetTitle("m_{\gamma e^{+}e^{-}} [GeV/c^{2}]");

    
    ROOT.SetOwnership(frame1,False);
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    # if option == 4 :
    #     leg.SetTextSize(0.02); 
    # if option == 5: 
    #     leg.SetTextSize(0.025); 
    # else:
    leg.SetTextSize(0.03);
    leg.Draw("");
    ROOT.SetOwnership(leg,False);    
    # if option == 4:
    #     txt = TPaveText(0.90,0.65,0.9,0.75,"NDC");
    if option == 3 or option == 5 or option == 4 or option == 25: 
        txt = TPaveText(0.90,0.85,0.9,0.95,"NDC");
    elif option == 30 : 
        txt = TPaveText(0.12,0.85,0.45,0.95,"NDC");
    elif option == 20: 
        txt = TPaveText(0.92,0.85,0.92,0.95,"NDC");
    elif option == 0 or option == 10 or option == 21:
        txt = TPaveText(0.1,0.85,0.29,0.95,"NDC");
    else:
        txt = TPaveText(0.2,0.85,0.45,0.95,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(33);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.03);
    txt.AddText("Simulation this thesis");
    txt.Draw();
    ROOT.SetOwnership(txt,False);
    # if option == 4:
    #     txt2 = TPaveText(0.845,0.6,0.83,0.725,"NDC");
    # if option == 3 or option == 5 or option == 4 or option == 25: 
    #     txt2 = TPaveText(0.845,0.8,0.83,0.925,"NDC");
    # elif option == 30: 
    #     txt2 = TPaveText(0.1,0.8,0.35,0.925,"NDC");
    # elif option == 20: 
    #     txt2 = TPaveText(0.89,0.8,0.89, 0.925,"NDC");
    # elif option == 0 or option == 10: 
    #     txt2 = TPaveText(0.1,0.8,0.27, 0.925,"NDC");
    # else:
    #     txt2 = TPaveText(0.2,0.8,0.35,0.925,"NDC");
    # txt2.SetFillColor(kWhite);
    # txt2.SetFillStyle(0);
    # txt2.SetBorderSize(0);
    # txt2.SetTextAlign(33);#middle,left
    # txt2.SetTextFont(42);#helvetica
    # txt2.SetTextSize(0.03);
    # txt2.AddText("this thesis");
    # txt2.Draw();
    # ROOT.SetOwnership(txt2,False);
    # # if option == 4:
    # #     txt3 = TPaveText(0.9,0.55,0.9,0.70,"NDC");
    # if option == 3 or option == 5 or option == 4 or option == 25: 
    #     txt3 = TPaveText(0.9,0.75,0.9,0.90,"NDC");
    # elif option == 30: 
    #     txt3 = TPaveText(0.1,0.75,0.4,0.90,"NDC");
    # elif option == 20 : 
    #     txt3 = TPaveText(0.92,0.75,0.92,0.90,"NDC");
    # elif option == 0 or option == 10: 
    #     txt3 = TPaveText(0.1,0.75,0.3,0.90,"NDC");
    # else:
    #     txt3 = TPaveText(0.2,0.75,0.4,0.90,"NDC");
    
    if option == 3 or option == 5 or option == 4 or option == 25: 
        txt3 = TPaveText(0.88,0.8,0.88,0.925,"NDC");
    elif option == 30 :
        txt3 = TPaveText(0.125,0.8,0.425,0.925,"NDC");
    elif option == 20: 
        txt3 = TPaveText(0.905,0.8,0.905, 0.925,"NDC");
    elif option == 0 or option == 10 or option == 21: 
        txt3 = TPaveText(0.1,0.8,0.275, 0.925,"NDC");
    else:
        txt3 = TPaveText(0.2,0.8,0.425,0.925,"NDC");
    txt3.SetFillColor(kWhite);
    txt3.SetFillStyle(0);
    txt3.SetBorderSize(0);
    txt3.SetTextAlign(33);#middle,left
    txt3.SetTextFont(42);#helvetica
    txt3.SetTextSize(0.03);
    txt3.AddText("pp, #sqrt{s} = 13.6TeV");
    txt3.Draw();
    ROOT.SetOwnership(txt3,False);

    if (option == 1 or option == 2 ):
        txt5 = TPaveText(0.2,0.75,0.4,0.9,"NDC");
        txt5.SetFillColor(kWhite);
        txt5.SetFillStyle(0);
        txt5.SetBorderSize(0);
        txt5.SetTextAlign(33);#middle,left
        txt5.SetTextFont(42);#helvetica
        txt5.SetTextSize(0.03);
        txt5.AddText("MC matched");
        txt5.Draw();
        ROOT.SetOwnership(txt5,False);
    if option == 0: 
        txt5 = TPaveText(0.15,0.75,0.25,0.9,"NDC");
        txt5.SetFillColor(kWhite);
        txt5.SetFillStyle(0);
        txt5.SetBorderSize(0);
        txt5.SetTextAlign(33);#middle,left
        txt5.SetTextFont(42);#helvetica
        txt5.SetTextSize(0.03);
        txt5.AddText("MC matched");
        txt5.Draw();
        ROOT.SetOwnership(txt5,False);
    if option == 10: 
        txt5 = TPaveText(0.1,0.75,0.26,0.9,"NDC");
        txt5.SetFillColor(kWhite);
        txt5.SetFillStyle(0);
        txt5.SetBorderSize(0);
        txt5.SetTextAlign(33);#middle,left
        txt5.SetTextFont(42);#helvetica
        txt5.SetTextSize(0.03);
        txt5.AddText("MC generated");
        txt5.Draw();
        ROOT.SetOwnership(txt5,False);
    if option == 5:
        txt5 = TPaveText(0.7,0.75,0.865,0.9,"NDC");
        txt5.SetFillColor(kWhite);
        txt5.SetFillStyle(0);
        txt5.SetBorderSize(0);
        txt5.SetTextAlign(33);#middle,left
        txt5.SetTextFont(42);#helvetica
        txt5.SetTextSize(0.03);
        txt5.AddText("MC matched");
        txt5.Draw();
        ROOT.SetOwnership(txt5,False);
    if (option == 30):
        txt4 = TPaveText(0.35,0.75,0.395,0.9,"NDC");
        txt4.SetFillColor(kWhite);
        txt4.SetFillStyle(0);
        txt4.SetBorderSize(0);
        txt4.SetTextAlign(33);#middle,left
        txt4.SetTextFont(42);#helvetica
        txt4.SetTextSize(0.03);
        txt4.AddText("MC generated");
        txt4.Draw();
        ROOT.SetOwnership(txt4,False);
    # txt6 = TPaveText(0.1,0.9,0.4,0.90,"NDC");
    # txt6.SetFillColor(kWhite);
    # txt6.SetFillStyle(0);
    # txt6.SetBorderSize(0);
    # txt6.SetTextAlign(33);#middle,left
    # txt6.SetTextFont(42);#helvetica
    # txt6.SetTextSize(0.03);
    # txt6.AddText("|\eta| < 0.9");
    # txt6.Draw();
    # ROOT.SetOwnership(txt6,False);


    #Draw arrows at the locations of J/psi, chic1 and chic2
    if option == 11:
        arrow = TArrow( 3.0969, 5000, 3.0969, ymin, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
    elif option == 10:
        arrow = TArrow( 3.0969, 250000, 3.0969, ymin, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
    elif option == 0:
        arrow = TArrow( 3.0969, 1000000, 3.0969, ymin, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
    elif option == 21:
        arrow = TArrow( 3.0969, 500000, 3.0969, ymin, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
    elif option == 25:
        arrow = TArrow( 3.0969, 0.8, 3.0969, ymin, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
    elif option == 20:
        arrow = TArrow( 3.0969, 1000000, 3.0969, ymin, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
    elif option == 22:
        arrow = TArrow( 3.0969, 1200000, 3.0969, ymin, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
    elif option == 30:
        arrow = TArrow( 3.51069, 2000000, 3.51069, ymin, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
        arrow2 = TArrow( 3.55617, 2000000, 3.55617, ymin, 0.02, '|>' )
        arrow2.SetFillStyle( 1001 )
        arrow2.Draw()
    elif option == 3:
        arrow = TArrow( 3.51069, ymin+50000, 3.51069, ymin, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
        arrow2 = TArrow( 3.55617, ymin+50000, 3.55617, ymin, 0.02, '|>' )
        arrow2.SetFillStyle( 1001 )
        arrow2.Draw()
    elif option == 5:
        arrow = TArrow( 3.51069, 1500, 3.51069, ymin, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
        arrow2 = TArrow( 3.55617, 1500, 3.55617, ymin, 0.02, '|>' )
        arrow2.SetFillStyle( 1001 )
        arrow2.Draw()
    elif option == 4: 
        arrow = TArrow( 3.51069, 1500, 3.51069, ymin, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
        arrow2 = TArrow( 3.55617, 1500, 3.55617, ymin, 0.02, '|>' )
        arrow2.SetFillStyle( 1001 )
        arrow2.Draw()

    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs(plotname);



if __name__ == "__main__":
    
    #mass_plot("AnalysisResults_chic_20240215.root", 11, "plot_mass_eejpsi_truthmatched.svg")
    #mass_plot("AnalysisResults_chic_20240217.root", 15, "20240217/plot_mass_ee_eejpsi_matched.pdf")´

    filename = "AnalysisResults_chicall_20240224.root"
    # mass_plot(filename, 11, "20240305/plot_mass_eejpsi_truthmatched.pdf") 
    # mass_plot(filename, 11, "20240305/plot_mass_eejpsi_truthmatched.svg")  
    mass_plot(filename, 10, "20240305/plot_mass_eejpsi_truth.pdf") 
    mass_plot(filename, 10, "20240305/plot_mass_eejpsi_truth.svg") 
    mass_plot(filename, 0, "20240305/plot_mass_eejpsi_matched.pdf") 
    mass_plot(filename, 0, "20240305/plot_mass_eejpsi_matched.svg") 
    mass_plot(filename, 30, "20240305/plot_mass_eephotonchic_truth.pdf") 
    mass_plot(filename, 30, "20240305/plot_mass_eephotonchic_truth.svg") 
    mass_plot(filename, 21, "20240305/plot_mass_eewithoutbackground_eejpsi_matched.pdf") 
    mass_plot(filename, 21, "20240305/plot_mass_eewithoutbackground_eejpsi_matched.svg") 
    mass_plot(filename, 3, "20240305/plot_mass_eephoton_triple.pdf") 
    mass_plot(filename, 3, "20240305/plot_mass_eephoton_triple.svg") 
    mass_plot(filename, 5, "20240305/plot_mass_eephotonchic_matched.pdf") 
    mass_plot(filename, 5, "20240305/plot_mass_eephotonchic_matched.svg") 
    mass_plot(filename, 4, "20240305/plot_mass_eephotonchic_triplematched.pdf") 
    mass_plot(filename, 4, "20240305/plot_mass_eephotonchic_triplematched.svg") 
    mass_plot(filename, 20, "20240305/plot_mass_ee_background.pdf") 
    mass_plot(filename, 20, "20240305/plot_mass_ee_background.svg") 
    mass_plot(filename, 22, "20240305/plot_mass_ee_background_smallrange.pdf") 
    mass_plot(filename, 22, "20240305/plot_mass_ee_background_smallrange.svg") 
    # mass_plot(filename, 25, "20240305/plot_ratiomass_triple_matched.pdf")
