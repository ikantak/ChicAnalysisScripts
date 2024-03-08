import re
import numpy as np
import datetime
import math

import ROOT
import ctypes
from ctypes import *
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency
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

def pT_plot(filename1, option, plotname):
    rootfile_data = TFile.Open(filename1, "READ");
    list_data = rootfile_data.Get("analysis-dilepton-photon");
    list_data2 = list_data.Get("output");

    if (option == 0 or option == 4 or option == 12 or option == 30):
        mc_cut_chic1 = list_data2.FindObject("MCTruthGen_cut_Chic1")
        mc_cut_chic2 = list_data2.FindObject("MCTruthGen_cut_Chic2") 

        chic1_cut_pt = mc_cut_chic1.FindObject("Pt");
        chic1_cut_pt.SetDirectory(0);
        make_common_style(chic1_cut_pt, kFullSquare, 1.0, kGray+3, 1, 0); #kViolet+2     #kRed      #kCyan+3
        ROOT.SetOwnership(chic1_cut_pt, False);
        print(chic1_cut_pt.GetNbinsX())
        chic1_cut_pt.Rebin(10)
        print(chic1_cut_pt.GetNbinsX())
        chic1_cut_pt.Scale(1, "width")

        chic2_cut_pt = mc_cut_chic2.FindObject("Pt");
        chic2_cut_pt.SetDirectory(0);
        make_common_style(chic2_cut_pt, kFullSquare, 1.0, kGray+1, 1, 0);#kMagenta+2    #kOrange+7     #kTeal-1
        ROOT.SetOwnership(chic2_cut_pt, False);
        print(chic2_cut_pt.GetNbinsX())
        chic2_cut_pt.Rebin(10)
        print(chic2_cut_pt.GetNbinsX())
        chic2_cut_pt.Scale(1, "width")
    
    if (option == 1 or option == 3 or option ==11): 
        mc_cut_jpsi = list_data2.FindObject("MCTruthGen_cut_Jpsi")
        jpsi_cut_pt = mc_cut_jpsi.FindObject("Pt");
        jpsi_cut_pt.SetDirectory(0);
        make_common_style(jpsi_cut_pt, kFullTriangleUp, 1.0, kBlue+2, 1, 0);#kSpring-7    #kGreen+3
        ROOT.SetOwnership(jpsi_cut_pt, False);
        print(jpsi_cut_pt.GetNbinsX())
        jpsi_cut_pt.Rebin(10)
        print(jpsi_cut_pt.GetNbinsX())
        jpsi_cut_pt.Scale(1, "width")
    
    if (option == 2 or option == 3 or option == 4): 
        mc_cut_jpsi_chic1 = list_data2.FindObject("MCTruthGen_cut_JpsiFromChic1")
        mc_cut_jpsi_chic2 = list_data2.FindObject("MCTruthGen_cut_JpsiFromChic2")
        jpsi1_cut_pt = mc_cut_jpsi_chic1.FindObject("Pt");
        jpsi1_cut_pt.SetDirectory(0);
        make_common_style(jpsi1_cut_pt, kFullTriangleUp, 1.0, kAzure-5, 1, 0);#kSpring-7    #kGreen+3
        ROOT.SetOwnership(jpsi1_cut_pt, False);
        print(jpsi1_cut_pt.GetNbinsX())
        jpsi1_cut_pt.Rebin(10)
        print(jpsi1_cut_pt.GetNbinsX())
        jpsi1_cut_pt.Scale(1, "width")

        jpsi2_cut_pt = mc_cut_jpsi_chic2.FindObject("Pt")
        jpsi2_cut_pt.SetDirectory(0);
        make_common_style(jpsi2_cut_pt, kFullTriangleUp, 1.0, kAzure+8, 1, 0);#kGreen+2  kSpring-8
        ROOT.SetOwnership(jpsi2_cut_pt, False);
        print(jpsi2_cut_pt.GetNbinsX())
        jpsi2_cut_pt.Rebin(10)
        print(jpsi2_cut_pt.GetNbinsX())
        jpsi2_cut_pt.Scale(1, "width")
    
    
    if (option == 5 or option == 25 or option == 11 or option ==23 or option == 6): 
        mc_cut_ee_jpsi = list_data2.FindObject("MCTruthGenPair_cut_eeFromJpsi")
        eejpsi_cut_pt = mc_cut_ee_jpsi.FindObject("Pt")
        eejpsi_cut_pt.SetDirectory(0)
        make_common_style(eejpsi_cut_pt, kFullCross, 1.0, kBlue+2, 1, 0)
        ROOT.SetOwnership(eejpsi_cut_pt, False)
        print(eejpsi_cut_pt.GetNbinsX())
        eejpsi_cut_pt.Rebin(10)
        print(eejpsi_cut_pt.GetNbinsX())
        eejpsi_cut_pt.Scale(1, "width")
        #eejpsi_cut_pt.Rebin(10)


    if (option == 6 or option ==26 or option == 12 or option ==23 or option == 30): 
        mc_cut_ee_chic1 = list_data2.FindObject("MCTruthGenPair_cut_eeFromJpsiFromChic1")
        mc_cut_ee_chic2 = list_data2.FindObject("MCTruthGenPair_cut_eeFromJpsiFromChic2")
        eechic1_cut_pt = mc_cut_ee_chic1.FindObject("Pt")
        eechic1_cut_pt.SetDirectory(0)
        make_common_style(eechic1_cut_pt, kFullCross, 1.0, kAzure-5, 1, 0)
        ROOT.SetOwnership(eechic1_cut_pt, False)
        

        eechic2_cut_pt = mc_cut_ee_chic2.FindObject("Pt")
        eechic2_cut_pt.SetDirectory(0)
        make_common_style(eechic2_cut_pt, kFullCross, 1.0, kAzure+8, 1, 0)
        ROOT.SetOwnership(eechic2_cut_pt, False)
        print(eechic1_cut_pt.GetNbinsX())
        eechic1_cut_pt.Rebin(10)
        print(eechic1_cut_pt.GetNbinsX())
        print(eechic2_cut_pt.GetNbinsX())
        eechic2_cut_pt.Rebin(10)
        print(eechic2_cut_pt.GetNbinsX())
        eechic1_cut_pt.Scale(1, "width")
        eechic2_cut_pt.Scale(1, "width")
        #eechic1_cut_pt.Rebin(10)
        #eechic2_cut_pt.Rebin(10)
    
    if (option == 7 or option == 27 or option ==12 or option == 30): 
        mc_cut_photon_chic1 = list_data2.FindObject("MCTruthGen_cut_PhotonFromChic1")
        mc_cut_photon_chic2 = list_data2.FindObject("MCTruthGen_cut_PhotonFromChic2")
        photon1_cut_pt = mc_cut_photon_chic1.FindObject("PtMC_photon");

        photon1_cut_pt.SetDirectory(0);
        make_common_style(photon1_cut_pt, kFullCircle, 1.0, kGreen+3, 1, 0);
        ROOT.SetOwnership(photon1_cut_pt, False);
        photon1_cut_pt.Rebin(50)
        photon1_cut_pt.Scale(1, "width")

        photon2_cut_pt = mc_cut_photon_chic2.FindObject("PtMC_photon");
        photon2_cut_pt.SetDirectory(0);
        make_common_style(photon2_cut_pt, kFullCircle, 1.0, kSpring-6, 1, 0);
        ROOT.SetOwnership(photon2_cut_pt, False);
        photon2_cut_pt.Rebin(50)
        photon2_cut_pt.Scale(1, "width")
        print(photon1_cut_pt.GetNbinsX())
        print(photon2_cut_pt.GetNbinsX())

    if (option == 8 or option == 28 or option == 30): 
        mc_cut_eephoton_chic1 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic1")
        mc_cut_eephoton_chic2 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic2")

        eephoton1_cut_pt = mc_cut_eephoton_chic1.FindObject("Pt_DileptonPhoton");
        eephoton1_cut_pt.SetDirectory(0);
        make_common_style(eephoton1_cut_pt, kFullCrossX, 1.0, kRed, 1, 0);
        ROOT.SetOwnership(eephoton1_cut_pt, False);
        eephoton1_cut_pt.Scale(1, "width")

        eephoton2_cut_pt = mc_cut_eephoton_chic2.FindObject("Pt_DileptonPhoton");
        eephoton2_cut_pt.SetDirectory(0);
        make_common_style(eephoton2_cut_pt, kFullCrossX, 1.0, kOrange+7, 1, 0);
        ROOT.SetOwnership(eephoton2_cut_pt, False);
        eephoton2_cut_pt.Scale(1, "width")

    if (option == 9 or option == 29): 
        mc_cut_eephoton_chic1 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic1")
        mc_cut_eephoton_chic2 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic2")
        mc_cut_eephoton_chic12 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic12")

        eephoton1_cut_pt = mc_cut_eephoton_chic1.FindObject("Pt_DileptonPhoton");
        eephoton1_cut_pt.SetDirectory(0);
        make_common_style(eephoton1_cut_pt, kFullCrossX, 1.0, kRed, 1, 0);
        ROOT.SetOwnership(eephoton1_cut_pt, False);
        eephoton12_cut_pt.Scale(1, "width")

        eephoton2_cut_pt = mc_cut_eephoton_chic2.FindObject("Pt_DileptonPhoton");
        eephoton2_cut_pt.SetDirectory(0);
        make_common_style(eephoton2_cut_pt, kFullCrossX, 1.0, kOrange+7, 1, 0);
        ROOT.SetOwnership(eephoton2_cut_pt, False);
        eephoton2_cut_pt.Scale(1, "width")

        eephoton12_cut_pt = mc_cut_eephoton_chic12.FindObject("Pt_DileptonPhoton");
        eephoton12_cut_pt.SetDirectory(0);
        make_common_style(eephoton12_cut_pt, kFullCrossX, 1.0, kRed+2, 1, 0);
        ROOT.SetOwnership(eephoton12_cut_pt, False);
        eephoton12_cut_pt.Scale(1, "width")
    
    
    if (option == 15 or option == 25 or option ==13 or option == 23):
        data_cut_ee_jpsi = list_data2.FindObject("DileptonsSelected_cut_matchedMC_eeFromJpsi")
        eejpsi_cut_data_pt = data_cut_ee_jpsi.FindObject("Pt")
        eejpsi_cut_data_pt.SetDirectory(0)
        make_common_style(eejpsi_cut_data_pt, kFullCross, 1.0, kCyan+3, 1, 0)
        ROOT.SetOwnership(eejpsi_cut_data_pt, False)
        print(eejpsi_cut_data_pt.GetNbinsX())
        eejpsi_cut_data_pt.Rebin(10)
        print(eejpsi_cut_data_pt.GetNbinsX())
        eejpsi_cut_data_pt.Scale(1, "width")    
    
    if (option == 16 or option == 26 or option == 13 or option == 23 or option == 22 or option == 30):
        data_cut_ee_chic1 = list_data2.FindObject("DileptonsSelected_cut_matchedMC_eeFromJpsiFromChic1")
        data_cut_ee_chic2 = list_data2.FindObject("DileptonsSelected_cut_matchedMC_eeFromJpsiFromChic2")

        eechic1_cut_data_pt = data_cut_ee_chic1.FindObject("Pt");
        eechic1_cut_data_pt.SetDirectory(0);
        make_common_style(eechic1_cut_data_pt, kFullCross, 1.0, kCyan+1, 1, 0);
        ROOT.SetOwnership(eechic1_cut_data_pt, False);
        
        eechic2_cut_data_pt = data_cut_ee_chic2.FindObject("Pt");
        eechic2_cut_data_pt.SetDirectory(0);
        make_common_style(eechic2_cut_data_pt, kFullCross, 1.0, kCyan-9, 1, 0);
        ROOT.SetOwnership(eechic2_cut_data_pt, False);
        print(eechic1_cut_data_pt.GetNbinsX())
        eechic1_cut_data_pt.Rebin(10)
        print(eechic1_cut_data_pt.GetNbinsX())
        print(eechic2_cut_data_pt.GetNbinsX())
        eechic2_cut_data_pt.Rebin(10)
        print(eechic2_cut_data_pt.GetNbinsX())
        eechic1_cut_data_pt.Scale(1, "width")
        eechic2_cut_data_pt.Scale(1, "width")
    
    if (option == 17 or option == 27 or option == 12 or option == 22 or option == 30):
        data_cut_photon_chic1 = list_data2.FindObject("Selected_cut_matchedMC_PhotonFromChic1")
        data_cut_photon_chic2 = list_data2.FindObject("Selected_cut_matchedMC_PhotonFromChic2")
        photonchic1_cut_data_pt = data_cut_photon_chic1.FindObject("Pt_Photon");
        print(photonchic1_cut_data_pt.GetNbinsX())
        photonchic1_cut_data_pt.Rebin(50)
        print(photonchic1_cut_data_pt.GetNbinsX())
        photonchic1_cut_data_pt.SetDirectory(0);
        make_common_style(photonchic1_cut_data_pt, kFullCircle, 1.0, kGreen+1, 1, 0);#kViolet+2     #kRed      #kCyan+3
        ROOT.SetOwnership(photonchic1_cut_data_pt, False);
        photonchic1_cut_data_pt.Scale(1, "width")

        photonchic2_cut_data_pt = data_cut_photon_chic2.FindObject("Pt_Photon");
        photonchic2_cut_data_pt.Rebin(50)
        photonchic2_cut_data_pt.SetDirectory(0);
        make_common_style(photonchic2_cut_data_pt, kFullCircle, 1.0, kSpring+7, 1, 0);#kMagenta+2    #kOrange+7     #kTeal-1
        ROOT.SetOwnership(photonchic2_cut_data_pt, False)
        photonchic2_cut_data_pt.Scale(1, "width")
        print(photonchic2_cut_data_pt.GetNbinsX())

    if (option == 18 or option == 19 or option == 28 or option == 29 or option == 22 or option == 30): 
        data_cut_eephoton_chic1 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic1")
        data_cut_eephoton_chic2 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic2")

        eephoton1_cut_data_pt = data_cut_eephoton_chic1.FindObject("Pt_DileptonPhoton");
        eephoton1_cut_data_pt.SetDirectory(0);
        make_common_style(eephoton1_cut_data_pt, kFullCrossX, 1.0, kPink+7, 1, 0);
        ROOT.SetOwnership(eephoton1_cut_data_pt, False);
        eephoton1_cut_data_pt.Rebin(100)
        eephoton1_cut_data_pt.Scale(1, "width")

        eephoton2_cut_data_pt = data_cut_eephoton_chic2.FindObject("Pt_DileptonPhoton");
        eephoton2_cut_data_pt.SetDirectory(0);
        make_common_style(eephoton2_cut_data_pt, kFullCrossX, 1.0, kPink+1, 1, 0);
        ROOT.SetOwnership(eephoton2_cut_data_pt, False)
        eephoton2_cut_data_pt.Rebin(100)
        eephoton2_cut_data_pt.Scale(1, "width")

    if (option == 19 or option == 29 ):
        data_cut_eephoton_chic12 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic12")
        eephoton12_cut_data_pt = data_cut_eephoton_chic12.FindObject("Pt_DileptonPhoton");
        eephoton12_cut_data_pt.SetDirectory(0);
        make_common_style(eephoton12_cut_data_pt, kFullCrossX, 1.0, kPink-8, 1, 0);
        ROOT.SetOwnership(eephoton12_cut_data_pt, False);
        eephoton12_cut_data_pt.Scale(1, "width")

    
    ymax = 0
    if option == 0 : 
        ymax= max(chic1_cut_pt.GetMaximum(), chic2_cut_pt.GetMaximum())
        ymin = min(chic1_cut_pt.GetMinimum(), chic2_cut_pt.GetMinimum())
    elif option == 1 : 
        ymax = jpsi_cut_pt.GetMaximum()
        ymin = jpsi_cut_pt.GetMinimum()
    elif option == 2: 
        ymax = max(jpsi1_cut_pt.GetMaximum(), jpsi2_cut_pt.GetMaximum())
        ymin = min(jpsi1_cut_pt.GetMinimum(), jpsi2_cut_pt.GetMinimum())
    elif option == 3: 
        ymax = max(jpsi1_cut_pt.GetMaximum(), jpsi2_cut_pt.GetMaximum(),jpsi_cut_pt.GetMaximum())
        ymin = min(jpsi1_cut_pt.GetMinimum(), jpsi2_cut_pt.GetMinimum(), jpsi_cut_pt.GetMinimum())
    elif option == 4: 
        ymax= max(chic1_cut_pt.GetMaximum(), chic2_cut_pt.GetMaximum(), jpsi1_cut_pt.GetMaximum(), jpsi2_cut_pt.GetMaximum())
        ymin = min(chic1_cut_pt.GetMinimum(), chic2_cut_pt.GetMinimum(), jpsi1_cut_pt.GetMinimum(), jpsi2_cut_pt.GetMinimum())
    elif option == 5: 
        ymax = eejpsi_cut_pt.GetMaximum()
        ymin = eejpsi_cut_pt.GetMinimum()
    elif option == 6: 
        ymax = max(eechic1_cut_pt.GetMaximum(), eechic2_cut_pt.GetMaximum(), eejpsi_cut_pt.GetMaximum())
        ymin = min(eechic1_cut_pt.GetMinimum(), eechic2_cut_pt.GetMinimum())
    elif option == 7: 
        ymax = max(photon1_cut_pt.GetMaximum(), photon2_cut_pt.GetMaximum())
        ymin = min(photon1_cut_pt.GetMinimum(), photon2_cut_pt.GetMinimum())
    elif option == 8: 
        ymax = max(eephoton1_cut_pt.GetMaximum(), eephoton2_cut_pt.GetMaximum())
        ymin = min(eephoton1_cut_pt.GetMinimum(), eephoton2_cut_pt.GetMinimum())
    elif option == 9: 
        ymax = max(eephoton1_cut_pt.GetMaximum(),  eephoton2_cut_pt.GetMaximum(),  eephoton12_cut_pt.GetMaximum())
        ymin = min(eephoton1_cut_pt.GetMinimum(), eephoton2_cut_pt.GetMinimum(), eephoton12_cut_pt.GetMinimum())
    
    elif option == 11:
        ymax = max(jpsi_cut_pt.GetMaximum(), eejpsi_cut_pt.GetMaximum())
        ymin = min(jpsi_cut_pt.GetMinimum(), eejpsi_cut_pt.GetMinimum())
    elif option == 12:
        ymax= max(chic1_cut_pt.GetMaximum(), chic2_cut_pt.GetMaximum(), photon1_cut_pt.GetMaximum(), photon2_cut_pt.GetMaximum())
        ymin = min(eechic1_cut_pt.GetMinimum(), eechic2_cut_pt.GetMinimum(), photon1_cut_pt.GetMinimum(), photon2_cut_pt.GetMinimum())
    elif option == 13:
        ymax = max(eejpsi_cut_data_pt.GetMaximum(), eechic1_cut_data_pt.GetMaximum(), eechic2_cut_data_pt.GetMaximum())*1.5
        ymin = min(eejpsi_cut_data_pt.GetMinimum(), eechic1_cut_data_pt.GetMinimum(), eechic2_cut_data_pt.GetMinimum())
    elif option == 15: 
        ymax = eejpsi_cut_data_pt.GetMaximum()
        ymin = eejpsi_cut_data_pt.GetMinimum()
    elif option == 16: 
        ymax = max(eechic1_cut_data_pt.GetMaximum(), eechic2_cut_data_pt.GetMaximum())
        ymin = min(eechic1_cut_data_pt.GetMinimum(), eechic2_cut_data_pt.GetMinimum())
    elif option == 17: 
        ymax = max(photonchic1_cut_data_pt.GetMaximum(), photonchic2_cut_data_pt.GetMaximum())
        ymin = min(photonchic1_cut_data_pt.GetMinimum(), photonchic2_cut_data_pt.GetMinimum())
    elif option == 18: 
        ymax = max(eephoton1_cut_data_pt.GetMaximum(), eephoton2_cut_data_pt.GetMaximum())*2
        ymin = min(eephoton1_cut_data_pt.GetMinimum(), eephoton2_cut_data_pt.GetMinimum())
    elif option == 19:
        ymax = max(eephoton1_cut_data_pt.GetMaximum(), eephoton2_cut_data_pt.GetMaximum(), eephoton12_cut_data_pt.GetMaximum())
        ymin = min(eephoton1_cut_data_pt.GetMinimum(), eephoton2_cut_data_pt.GetMinimum(), eephoton12_cut_data_pt.GetMinimum())
    elif option == 22:
        ymax = max( eechic1_cut_data_pt.GetMaximum(), eechic2_cut_data_pt.GetMaximum(), photonchic1_cut_data_pt.GetMaximum(), photonchic2_cut_data_pt.GetMaximum())
        ymin = min( eechic1_cut_data_pt.GetMinimum(), eechic2_cut_data_pt.GetMinimum(), photonchic1_cut_data_pt.GetMinimum(), photonchic2_cut_data_pt.GetMinimum())
    elif option == 23:
        ymax = max(eejpsi_cut_pt.GetMaximum(), eejpsi_cut_data_pt.GetMaximum(), eechic1_cut_pt.GetMaximum(), eechic2_cut_pt.GetMaximum(), eechic1_cut_data_pt.GetMaximum(), eechic2_cut_data_pt.GetMaximum())
        ymin = min(eejpsi_cut_pt.GetMinimum(), eejpsi_cut_data_pt.GetMinimum(), eechic1_cut_pt.GetMinimum(), eechic2_cut_pt.GetMinimum(), eechic1_cut_data_pt.GetMinimum(), eechic2_cut_data_pt.GetMinimum())
    elif option == 25:
        ymax = max(eejpsi_cut_pt.GetMaximum(), eejpsi_cut_data_pt.GetMaximum())
        ymin = min(eejpsi_cut_pt.GetMinimum(), eejpsi_cut_data_pt.GetMinimum())
    elif option == 26: 
        ymax = max(eechic1_cut_pt.GetMaximum(), eechic2_cut_pt.GetMaximum(), eechic1_cut_data_pt.GetMaximum(), eechic2_cut_data_pt.GetMaximum())
        ymin = max(eechic1_cut_pt.GetMinimum(), eechic2_cut_pt.GetMinimum(), eechic1_cut_data_pt.GetMinimum(), eechic2_cut_data_pt.GetMinimum())
    elif option == 27:
        ymax = max(photon1_cut_pt.GetMaximum(), photon2_cut_pt.GetMaximum(), photonchic1_cut_data_pt.GetMaximum(), photonchic2_cut_data_pt.GetMaximum())*1.5
        ymin = 10
    elif option == 28:
        ymax = max(eephoton1_cut_pt.GetMaximum(), eephoton2_cut_pt.GetMaximum(), eephoton1_cut_data_pt.GetMaximum(), eephoton2_cut_data_pt.GetMaximum())
        ymin = min(eephoton1_cut_pt.GetMinimum(), eephoton2_cut_pt.GetMinimum(), eephoton1_cut_data_pt.GetMinimum(), eephoton2_cut_data_pt.GetMinimum())
    elif option == 29: 
        ymax = max(eephoton1_cut_pt.GetMaximum(),  eephoton2_cut_pt.GetMaximum(),  eephoton12_cut_pt.GetMaximum())
        ymin = min(eephoton1_cut_data_pt.GetMinimum(), eephoton2_cut_data_pt.GetMinimum(), eephoton12_cut_data_pt.GetMinimum())
    elif option == 30: 
        ymax = max(eechic1_cut_pt.GetMaximum(), chic1_cut_pt.GetMaximum())
        ymin = eephoton1_cut_data_pt.GetMinimum()
    
    
    #ymax = max(jpsi1_pt.GetMaximum(), jpsi2_pt.GetMaximum(), photon1_pt.GetMaximum(), photon2_pt.GetMaximum());
    if option == 18: 
        c1 = TCanvas("pT_distribution","pT distribution",0,0,1500,900);
    else: 
        c1 = TCanvas("pT_distribution","pT distribution",0,0,900,900);
    p1 = c1.cd();
    p1.SetPad(0,0.01,1,1);
    if option == 18: 
        p1.SetMargin(0.1,0.05,0.12,0.05);
    else: 
        p1.SetMargin(0.15,0.05,0.12,0.05);
    p1.SetTicks(1,1);
   
    p1.SetLogy() #log scale
    if ymax <10:
        ymax += 2
    if (ymax >9 and ymax <100):
        ymax += 20
    if (ymax >=100 and ymax < 1000):
        ymax += 100
    if (ymax >= 1000 and ymax < 10000):
        ymax += 2000
    if (ymax >= 10000 and ymax < 100000):
        ymax += 20000
    if (ymax >= 100000 and ymax < 1000000):
        ymax += 200000
    if (ymax >= 1000000 and ymax < 10000000):
        ymax += 2000000

    if (ymin < 1 or (ymin >= 1 and ymin < 10)):
        ymin = 1
    elif (ymin >= 10 and ymin < 100):
        ymin = 50
    elif (ymin >= 100 and ymin < 1000):
        ymin = 200
    elif (ymin>= 1000 and ymin < 10000):
        ymin = 3000
    elif (ymin >= 10000 and ymin < 100000):
        ymin = 3000
    elif (ymin >= 100000 and ymin < 1000000):
        ymin = 50000
    else: 
        ymin = 500000
    
    
    xmax = 15
    if ((option >=17 and option<19) or (option >=27 and option <=29)):
        xmax = 5
    if option == 27:
        xmax = 4.5
    if option == 13 or option == 6 or option == 0 or option == 3: 
        xmax = 20
    if option == 18: 
        xmax = 20
       

    frame1 = p1.DrawFrame(0, ymin, xmax, ymax);
    frame1.GetXaxis().SetTitle("p_{T} [GeV/c]");
    frame1.GetYaxis().SetTitle("Counts per GeV/c");
    frame1.GetXaxis().SetTitleSize(0.045);
    frame1.GetYaxis().SetTitleSize(0.045);
    if option == 18: 
        frame1.GetXaxis().SetTitleOffset(1.1);
        frame1.GetYaxis().SetTitleOffset(1.);
    else: 
        frame1.GetXaxis().SetTitleOffset(1.1);
        frame1.GetYaxis().SetTitleOffset(1.4);
    frame1.GetXaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetLabelSize(0.035);
    frame1.GetYaxis().SetMaxDigits(3);
    frame1.GetXaxis().SetLabelOffset(0.01);
    frame1.GetYaxis().SetLabelOffset(0.01);
    ROOT.SetOwnership(frame1,False);

    #h1data.SetMarkerStyle(kFullCircle)
    if option == 0:
        leg = TLegend(0.2,0.2,0.5,0.3);
        chic1_cut_pt.Draw("Esame")
        chic2_cut_pt.Draw("Esame")
        leg.AddEntry(chic1_cut_pt, "\chi_{c1} ","LP");
        leg.AddEntry(chic2_cut_pt, "\chi_{c2} ","LP");

    if option == 1: 
        leg = TLegend(0.2,0.2,0.5,0.3);
        jpsi_cut_pt.Draw("Esame")
        leg.AddEntry(jpsi_cut_pt, "J/\psi ", "LP")

    if option == 2: 
        leg = TLegend(0.2,0.2,0.5,0.3);
        jpsi2_cut_pt.Draw("Esame")
        jpsi1_cut_pt.Draw("Esame")
        leg.AddEntry(jpsi1_cut_pt, "J/\psi from \chi_{c1} ","LP");
        leg.AddEntry(jpsi2_cut_pt, "J/\psi from \chi_{c2} ","LP");
    
    if option == 3:
        leg = TLegend(0.2,0.2,0.5,0.35);
        jpsi_cut_pt.Draw("Esame,hist")
        jpsi2_cut_pt.Draw("Esame")
        jpsi1_cut_pt.Draw("Esame")
        leg.AddEntry(jpsi_cut_pt, "J/\psi from \chi_{c1} and \chi_{c2} ", "LP")
        leg.AddEntry(jpsi1_cut_pt, "J/\psi from \chi_{c1} ","LP");
        leg.AddEntry(jpsi2_cut_pt, "J/\psi from \chi_{c2} ","LP");
    
    if option == 4: 
        leg = TLegend(0.2,0.2,0.5,0.3);
        chic1_cut_pt.Draw("Esame")
        chic2_cut_pt.Draw("Esame")
        jpsi2_cut_pt.Draw("Esame")
        jpsi1_cut_pt.Draw("Esame")
        leg.AddEntry(chic1_cut_pt, "\chi_{c1} ","LP");
        leg.AddEntry(chic2_cut_pt, "\chi_{c2} ","LP");
        leg.AddEntry(jpsi1_cut_pt, "J/\psi from \chi_{c1} ","LP");
        leg.AddEntry(jpsi2_cut_pt, "J/\psi from \chi_{c2} ","LP");
    
    if option == 5: 
        leg = TLegend(0.5,0.6,1.0,0.75);
        eejpsi_cut_pt.Draw("Esame")
        leg.AddEntry(eejpsi_cut_pt, "e^{+} e^{-} from J/\psi", "LP")

    if option == 6:
        leg = TLegend(0.65,0.6,1.0,0.7);
        eejpsi_cut_pt.Draw("Esame")
        eechic1_cut_pt.Draw("Esame")
        eechic2_cut_pt.Draw("Esame")
        leg.AddEntry(eejpsi_cut_pt, "J/\psi \\rightarrow e^{+} e^{-}", "LP")
        leg.AddEntry(eechic1_cut_pt, "e^{+} e^{-} from \chi_{c1} ", "LP")
        leg.AddEntry(eechic2_cut_pt, "e^{+} e^{-} from \chi_{c2} ", "LP")

    if option == 7:
        leg = TLegend(0.55,0.65,1.0,0.75);
        photon1_cut_pt.Draw("Esame")
        photon2_cut_pt.Draw("Esame")
        leg.AddEntry(photon1_cut_pt, "\gamma from \chi_{c1} ","LP");
        leg.AddEntry(photon2_cut_pt, "\gamma from \chi_{c2} ","LP");
    
    if option == 8: 
        leg = TLegend(0.5,0.6,1.0,0.75);
        eephoton1_cut_pt.Draw("Esame")
        eephoton2_cut_pt.Draw("Esame")
        leg.AddEntry(eephoton1_cut_pt, "\gamma e^{+} e^{-} from \chi_{c1}","LP");
        leg.AddEntry(eephoton2_cut_pt, "\gamma e^{+} e^{-} from \chi_{c2}","LP");
    
    if option == 9: 
        leg = TLegend(0.6,0.62,1.0,0.72);
        eephoton12_cut_pt.Draw("Esame")
        eephoton1_cut_pt.Draw("Esame")
        eephoton2_cut_pt.Draw("Esame")
        leg.AddEntry(eephoton12_cut_pt, "\gamma e^{+} e^{-} from \chi_{c1} and \chi_{c2}", "LP")
        leg.AddEntry(eephoton1_cut_pt, "\gamma e^{+} e^{-} from \chi_{c1} ","LP")
        leg.AddEntry(eephoton2_cut_pt, "\gamma e^{+} e^{-} from \chi_{c2}","LP")
        
    
    if option == 11:
        leg = TLegend(0.2,0.2,0.5,0.3);
        jpsi_cut_pt.Draw("Esame")
        leg.AddEntry(jpsi_cut_pt, "J/\psi MC truth", "LP")
        eejpsi_cut_pt.Draw("Esame")
        leg.AddEntry(eejpsi_cut_pt, "e^{+} e^{-} from J/\psi MC truth", "LP")

    if option == 12:
        leg = TLegend(0.15,0.2,0.45,0.4);
        chic1_cut_pt.Draw("Esame")
        chic2_cut_pt.Draw("Esame")
        eechic1_cut_pt.Draw("Esame")
        eechic2_cut_pt.Draw("Esame")
        photon1_cut_pt.Draw("Esame")
        photon2_cut_pt.Draw("Esame")
        leg.AddEntry(chic1_cut_pt, "\chi_{c1} MC truth","LP");
        leg.AddEntry(chic2_cut_pt, "\chi_{c2} MC truth","LP");
        leg.AddEntry(eechic1_cut_pt, "e^{+} e^{-} from \chi_{c1} MC generated", "LP")
        leg.AddEntry(eechic2_cut_pt, "e^{+} e^{-} from \chi_{c2} MC generated", "LP")
        leg.AddEntry(photon1_cut_pt, "\gamma from \chi_{c1} MC generated","LP");
        leg.AddEntry(photon2_cut_pt, "\gamma from \chi_{c2} MC generated","LP");
        
    if option == 13: 
        leg = TLegend(0.2,0.2,0.5,0.3);
        eejpsi_cut_data_pt.Draw("Esame")
        eechic1_cut_data_pt.Draw("Esame")
        eechic2_cut_data_pt.Draw("Esame")
        leg.AddEntry(eejpsi_cut_data_pt, "J/\psi \\rightarrow e^{+} e^{-}", "LP")
        leg.AddEntry(eechic1_cut_data_pt, "e^{+} e^{-} from \chi_{c1}", "LP")
        leg.AddEntry(eechic2_cut_data_pt, "e^{+} e^{-} from \chi_{c2}", "LP")

    if option == 15: 
        leg = TLegend(0.5,0.6,1.0,0.65);
        eejpsi_cut_data_pt.Draw("Esame")
        leg.AddEntry(eejpsi_cut_data_pt, "e^{+} e^{-} from J/\psi MC matched", "LP")

    if option == 16:
        leg = TLegend(0.5,0.6,1.0,0.7);
        eechic1_cut_data_pt.Draw("Esame")
        eechic2_cut_data_pt.Draw("Esame")
        leg.AddEntry(eechic1_cut_data_pt, "e^{+} e^{-} from \chi_{c1} MC matched", "LP")
        leg.AddEntry(eechic2_cut_data_pt, "e^{+} e^{-} from \chi_{c2} MC matched ", "LP")

    if option == 17: 
        leg = TLegend(0.55,0.6,1.0,0.7);
        photonchic1_cut_data_pt.Draw("Esame")
        photonchic2_cut_data_pt.Draw("Esame")
        leg.AddEntry(photonchic1_cut_data_pt, "\gamma from \chi_{c1} MC matched ","LP")
        leg.AddEntry(photonchic2_cut_data_pt, "\gamma from \chi_{c2} MC matched ","LP")

    if option == 18:
        leg = TLegend(0.13,0.2,0.35,0.3);
        eephoton1_cut_data_pt.Draw("Esame")
        eephoton2_cut_data_pt.Draw("Esame")
        leg.AddEntry(eephoton1_cut_data_pt, "\chi_{c1} \\rightarrow \gamma e^{+} e^{-}","LP")
        leg.AddEntry(eephoton2_cut_data_pt, "\chi_{c2} \\rightarrow \gamma e^{+} e^{-}","LP")

    if option == 19:
        leg = TLegend(0.35,0.65,1.0,0.75);
        eephoton12_cut_data_pt.Draw("Esame")
        eephoton1_cut_data_pt.Draw("Esame")
        eephoton2_cut_data_pt.Draw("Esame")
        leg.AddEntry(eephoton12_cut_data_pt, "\gamma e^{+} e^{-} from \chi_{c1} and \chi_{c2} matched MC","LP")
        leg.AddEntry(eephoton1_cut_data_pt, "\gamma e^{+} e^{-} from \chi_{c1} matched MC","LP")
        leg.AddEntry(eephoton2_cut_data_pt, "\gamma e^{+} e^{-} from \chi_{c2} matched MC","LP")
        
    if option == 23:
        leg = TLegend(0.5,0.55,1.0,0.75);
        eejpsi_cut_pt.Draw("Esame")
        eechic1_cut_pt.Draw("Esame")
        eechic2_cut_pt.Draw("Esame")
        leg.AddEntry(eejpsi_cut_pt, "e^{+} e^{-} from J/\psi, MC truth", "LP")
        leg.AddEntry(eechic1_cut_pt, "e^{+} e^{-} from \chi_{c1} MC truth", "LP")
        leg.AddEntry(eechic2_cut_pt, "e^{+} e^{-} from \chi_{c2} MC truth", "LP")
        eejpsi_cut_data_pt.Draw("Esame")
        eechic1_cut_data_pt.Draw("Esame")
        eechic2_cut_data_pt.Draw("Esame")
        leg.AddEntry(eejpsi_cut_data_pt, "e^{+} e^{-} from J/\psi matched MC", "LP")
        leg.AddEntry(eechic1_cut_data_pt, "e^{+} e^{-} from \chi_{c1} matched MC", "LP")
        leg.AddEntry(eechic2_cut_data_pt, "e^{+} e^{-} from \chi_{c2} matched MC", "LP")

    if option == 25:
        leg = TLegend(0.5,0.6,1.0,0.7);
        eejpsi_cut_pt.Draw("Esame")
        eejpsi_cut_data_pt.Draw("Esame")
        leg.AddEntry(eejpsi_cut_pt, "e^{+} e^{-} from J/\psi, MC generated", "LP")
        leg.AddEntry(eejpsi_cut_data_pt, "e^{+} e^{-} from J/\psi MC matched", "LP")

    if option == 26:
        leg = TLegend(0.5,0.6,1.0,0.75);
        eechic1_cut_pt.Draw("Esame")
        eechic2_cut_pt.Draw("Esame")
        leg.AddEntry(eechic1_cut_pt, "e^{+} e^{-} from \chi_{c1} MC truth", "LP")
        leg.AddEntry(eechic2_cut_pt, "e^{+} e^{-} from \chi_{c2} MC truth", "LP")
        eechic1_cut_data_pt.Draw("Esame")
        eechic2_cut_data_pt.Draw("Esame")
        leg.AddEntry(eechic1_cut_data_pt, "e^{+} e^{-} from \chi_{c1} matched MC", "LP")
        leg.AddEntry(eechic2_cut_data_pt, "e^{+} e^{-} from \chi_{c2} matched MC", "LP")

    if option == 27:
        leg = TLegend(0.53,0.67,0.93,0.82);
        photon1_cut_pt.Draw("Esame")
        photon2_cut_pt.Draw("Esame")
        photonchic1_cut_data_pt.Draw("Esame")
        photonchic2_cut_data_pt.Draw("Esame")
        leg.AddEntry(photon1_cut_pt, "\gamma from \chi_{c1} MC generated","LP");
        leg.AddEntry(photon2_cut_pt, "\gamma from \chi_{c2} MC generated","LP");
        leg.AddEntry(photonchic1_cut_data_pt, "\gamma from \chi_{c1} MC matched","LP")
        leg.AddEntry(photonchic2_cut_data_pt, "\gamma from \chi_{c2} MC matched","LP")


    if option == 28:
        leg = TLegend(0.5,0.6,1.0,0.75);
        eephoton1_cut_pt.Draw("Esame")
        eephoton2_cut_pt.Draw("Esame")
        eephoton1_cut_data_pt.Draw("Esame")
        eephoton2_cut_data_pt.Draw("Esame")      
        leg.AddEntry(eephoton1_cut_pt, "\gamma e^{+} e^{-} from \chi_{c1} MC truth","LP");
        leg.AddEntry(eephoton2_cut_pt, "\gamma e^{+} e^{-} from \chi_{c2} MC truth","LP");
        leg.AddEntry(eephoton1_cut_data_pt, "\gamma e^{+} e^{-} from \chi_{c1} matched MC","LP")
        leg.AddEntry(eephoton2_cut_data_pt, "\gamma e^{+} e^{-} from \chi_{c2} matched MC","LP")
    
    if option == 29: 
        leg = TLegend(0.2,0.25,0.5,0.45);
        eephoton12_cut_pt.Draw("Esame")
        eephoton1_cut_pt.Draw("Esame")
        eephoton2_cut_pt.Draw("Esame")
        leg.AddEntry(eephoton12_cut_pt, "\gamma e^{+} e^{-} from \chi_{c1} and \chi_{c2} MC truth", "LP")
        leg.AddEntry(eephoton1_cut_pt, "\gamma e^{+} e^{-} from \chi_{c1} MC truth","LP")
        leg.AddEntry(eephoton2_cut_pt, "\gamma e^{+} e^{-} from \chi_{c2} MC truth","LP")
        eephoton12_cut_data_pt.Draw("Esame")
        eephoton1_cut_data_pt.Draw("Esame")
        eephoton2_cut_data_pt.Draw("Esame")
        leg.AddEntry(eephoton12_cut_data_pt, "\gamma e^{+} e^{-} from \chi_{c1} and \chi_{c2} matched MC","LP")
        leg.AddEntry(eephoton1_cut_data_pt, "\gamma e^{+} e^{-} from \chi_{c1} matched MC","LP")
        leg.AddEntry(eephoton2_cut_data_pt, "\gamma e^{+} e^{-} from \chi_{c2} matched MC","LP")

    if option == 30: 
        leg = TLegend(0.2,0.25,0.5,0.45);
        chic1_cut_pt.Draw("Esame")
        leg.AddEntry(chic1_cut_pt, "\chi_{c1} MC truth","LP");
        eechic1_cut_pt.Draw("Esame")
        photon1_cut_pt.Draw("Esame")
        eephoton1_cut_pt.Draw("Esame")
        leg.AddEntry(eechic1_cut_pt, "e^{+} e^{-} from \chi_{c1} MC truth", "LP")
        leg.AddEntry(photon1_cut_pt, "\gamma from \chi_{c1} MC truth","LP");
        leg.AddEntry(eephoton1_cut_pt, "\gamma e^{+} e^{-} from \chi_{c1} MC truth","LP")
        # eechic1_data_pt.Draw("Esame")
        # photonchic1_data_pt.Draw("Esame")
        # eephoton1_data_pt.Draw("Esame")
        # leg.AddEntry(eechic1_data_pt, "e^{+} e^{-} from \chi_{c1} matched MC", "LP")
        # leg.AddEntry(photonchic1_data_pt, "\gamma from \chi_{c1} matched MC","LP")
        # leg.AddEntry(eephoton1_data_pt, "e^{+} e^{-} \gamma from \chi_{c1} matched MC","LP")

    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.03);
    leg.Draw("");
    ROOT.SetOwnership(leg,False);    
    if option == 18: 
        txt = TPaveText(0.92,0.85,0.92,0.95,"NDC");
    else: 
        txt = TPaveText(0.90,0.85,0.9,0.95,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(33);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.03);
    txt.AddText("Simulation this thesis");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    # if option == 18: 
    #     txt2 = TPaveText(0.89,0.8,0.89, 0.925,"NDC");
    # else: 
    #     txt2 = TPaveText(0.845,0.8,0.83,0.925,"NDC");
    # txt2.SetFillColor(kWhite);
    # txt2.SetFillStyle(0);
    # txt2.SetBorderSize(0);
    # txt2.SetTextAlign(33);#middle,left
    # txt2.SetTextFont(42);#helvetica
    # txt2.SetTextSize(0.03);
    # txt2.AddText("this thesis");
    # txt2.Draw();
    # ROOT.SetOwnership(txt2,False);

    # if option == 18: 
    #     txt3 = TPaveText(0.92,0.75,0.92,0.90,"NDC");
    # else: 
    #     txt3 = TPaveText(0.9,0.75,0.9,0.90,"NDC");
    if option == 18: 
        txt3 = TPaveText(0.905,0.8,0.905, 0.925,"NDC");
    else: 
        txt3 = TPaveText(0.88,0.8,0.88,0.925,"NDC");
    txt3.SetFillColor(kWhite);
    txt3.SetFillStyle(0);
    txt3.SetBorderSize(0);
    txt3.SetTextAlign(33);#middle,left
    txt3.SetTextFont(42);#helvetica
    txt3.SetTextSize(0.03);
    txt3.AddText("pp, #sqrt{s} = 13.6TeV");
    txt3.Draw();
    ROOT.SetOwnership(txt3,False);

    if (option <10):
        txt4 = TPaveText(0.84,0.75,0.863,0.9,"NDC");
        txt4.SetFillColor(kWhite);
        txt4.SetFillStyle(0);
        txt4.SetBorderSize(0);
        txt4.SetTextAlign(33);#middle,left
        txt4.SetTextFont(42);#helvetica
        txt4.SetTextSize(0.03);
        txt4.AddText("MC generated");
        txt4.Draw();
        ROOT.SetOwnership(txt4,False);
    elif (option >12 and option <23 and option != 18):
        txt5 = TPaveText(0.84,0.75,0.863,0.9,"NDC");
        txt5.SetFillColor(kWhite);
        txt5.SetFillStyle(0);
        txt5.SetBorderSize(0);
        txt5.SetTextAlign(33);#middle,left
        txt5.SetTextFont(42);#helvetica
        txt5.SetTextSize(0.03);
        txt5.AddText("MC matched");
        txt5.Draw();
        ROOT.SetOwnership(txt5,False);
    if option == 18: 
        txt5 = TPaveText(0.89,0.75,0.89,0.9,"NDC");
        txt5.SetFillColor(kWhite);
        txt5.SetFillStyle(0);
        txt5.SetBorderSize(0);
        txt5.SetTextAlign(33);#middle,left
        txt5.SetTextFont(42);#helvetica
        txt5.SetTextSize(0.03);
        txt5.AddText("MC matched");
        txt5.Draw();
        ROOT.SetOwnership(txt5,False);
    
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

    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs(plotname);



if __name__ == "__main__":
    filename = "AnalysisResults_chicall_20240224.root"
    pT_plot(filename, 0, "20240305/plot_etapt_chic12_truth.pdf")
    pT_plot(filename, 0, "20240305/plot_etapt_chic12_truth.svg")
    pT_plot(filename, 3, "20240305/plot_etapt_jpsichic12_truth.pdf")
    pT_plot(filename, 3, "20240305/plot_etapt_jpsichic12_truth.svg")
    pT_plot(filename, 6, "20240305/plot_etapt_eejpsichic12_truth.pdf")
    pT_plot(filename, 6, "20240305/plot_etapt_eejpsichic12_truth.svg")
    pT_plot(filename, 13, "20240305/plot_etapt_eejpsichic12_matched.pdf")
    pT_plot(filename, 13, "20240305/plot_etapt_eejpsichic12_matched.svg")
    pT_plot(filename, 27, "20240305/plot_etapt_photonchic12_truthmatched.pdf")
    pT_plot(filename, 27, "20240305/plot_etapt_photonchic12_truthmatched.svg")
    pT_plot(filename, 18, "20240305/plot_etapt_eephotonchic12_matched.pdf")
    pT_plot(filename, 18, "20240305/plot_etapt_eephotonchic12_matched.svg")

    
