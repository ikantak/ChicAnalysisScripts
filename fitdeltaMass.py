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

    if option == 0: 
        mc_eephoton_chic1 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic1")
        mc_eephoton_chic2 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic2")
        mc_eephoton_chic12 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic12")
        eephoton1_deltaM = mc_eephoton_chic1.FindObject("DeltaMass_Jpsi");
        eephoton1_deltaM.SetDirectory(0);
        make_common_style(eephoton1_deltaM, kFullCrossX, 1.0, kRed, 1, 0);
        ROOT.SetOwnership(eephoton1_deltaM, False);
        eephoton1_deltaM.Rebin(2)
        print(eephoton1_deltaM.GetNbinsX())
        eephoton1_deltaM.Scale(1, "width")

        eephoton2_deltaM = mc_eephoton_chic2.FindObject("DeltaMass_Jpsi");
        eephoton2_deltaM.SetDirectory(0);
        make_common_style(eephoton2_deltaM, kFullCrossX, 1.0, kOrange+7, 1, 0);
        ROOT.SetOwnership(eephoton2_deltaM, False);
        eephoton2_deltaM.Rebin(2)
        print(eephoton2_deltaM.GetNbinsX())
        eephoton2_deltaM.Scale(1, "width")
    
        eephoton12_deltaM = mc_eephoton_chic12.FindObject("DeltaMass_Jpsi");
        eephoton12_deltaM.SetDirectory(0);
        make_common_style(eephoton12_deltaM, kFullCrossX, 1.0, kRed+2, 1, 0);
        ROOT.SetOwnership(eephoton12_deltaM, False);
        eephoton12_deltaM.Rebin(2)
        print(eephoton12_deltaM.GetNbinsX())
        eephoton12_deltaM.Scale(1, "width")
    
    if option == 1 or option == 4 or option == 5 or option == 6: 
        data_eephoton_chic1 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic1")
        data_eephoton_chic2 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic2")
        eephoton1_data_deltaM = data_eephoton_chic1.FindObject("DeltaMass_Jpsi");
        eephoton1_data_deltaM.SetDirectory(0);
        make_common_style(eephoton1_data_deltaM, kFullCrossX, 1.0, kViolet-1, 1, 0);
        ROOT.SetOwnership(eephoton1_data_deltaM, False);
        eephoton1_data_deltaM.Rebin(5)
        eephoton1_data_deltaM.Scale(1, "width")

        eephoton2_data_deltaM = data_eephoton_chic2.FindObject("DeltaMass_Jpsi");
        eephoton2_data_deltaM.SetDirectory(0);
        make_common_style(eephoton2_data_deltaM, kFullCrossX, 1.0, kViolet-4, 1, 0);
        ROOT.SetOwnership(eephoton2_data_deltaM, False);
        eephoton2_data_deltaM.Rebin(5)
        print(eephoton2_data_deltaM.GetNbinsX())
        eephoton2_data_deltaM.Scale(1, "width")
    
        data_eephoton_chic12 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic12")
        eephoton12_data_deltaM = data_eephoton_chic12.FindObject("DeltaMass_Jpsi");
        eephoton12_data_deltaM.SetDirectory(0);
        make_common_style(eephoton12_data_deltaM, kFullCrossX, 1.0, kMagenta+3, 1, 0);
        ROOT.SetOwnership(eephoton12_data_deltaM, False);
        eephoton12_data_deltaM.Rebin(5)
        eephoton12_data_deltaM.Scale(1, "width")

    if option == 2: 
        mc_eephoton_chic1 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic1")
        eephoton1_deltaM = mc_eephoton_chic1.FindObject("DeltaMass_Jpsi");
        eephoton1_deltaM.SetDirectory(0);
        make_common_style(eephoton1_deltaM, kFullCrossX, 1.0, kRed, 1, 0);
        ROOT.SetOwnership(eephoton1_deltaM, False);

    if option == 3: 
        mc_eephoton_chic2 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic2")
        eephoton2_deltaM = mc_eephoton_chic2.FindObject("DeltaMass_Jpsi");
        eephoton2_deltaM.SetDirectory(0);
        make_common_style(eephoton2_deltaM, kFullCrossX, 1.0, kOrange+7, 1, 0);
        ROOT.SetOwnership(eephoton2_deltaM, False);
    
    if (option == 2 or option == 3):
        g2 = ROOT.TF1("g2", "gaus", 3.4, 3.6)
        g3 = ROOT.TF1("g3", "gaus", 3.5, 3.6)
        g2.SetParLimits(1, 3.4, 3.6)
        g3.SetParLimits(1, 3.5, 3.7)
        total = ROOT.TF1("total", "gaus(0)+gaus(3)", 3.4, 3.7)
    elif option == 1: 
        g2 = ROOT.TF1("g2", "gaus", 3.485, 3.52)
        g3 = ROOT.TF1("g3", "gaus", 3.53, 3.57)
        g2.SetParLimits(1, 3.485, 3.52)
        g3.SetParLimits(1, 3.53, 3.57)
        total = ROOT.TF1("total", "gaus(0)+gaus(3)", 3.48, 3.57)
    elif option==0:
        # g2 = ROOT.TF1("g2", "gaus", 3.4, 3.6)
        # g3 = ROOT.TF1("g3", "gaus", 3.52, 3.6)
        # g2.SetParLimits(1, 3.4, 3.52)
        # g3.SetParLimits(1, 3.52, 3.7)
        # total = ROOT.TF1("total", "gaus(0)+gaus(3)", 3.4, 3.7)
        g2 = ROOT.TF1("lorentz1", "(0.5*[0]*[1]/TMath::Pi()) / TMath::Max(1.e-10,(x-[2])*(x-[2])+ .25*[1]*[1])", 3.4, 3.6)
        g3 = ROOT.TF1("lorentz2", "(0.5*[0]*[1]/TMath::Pi()) / TMath::Max(1.e-10,(x-[2])*(x-[2])+ .25*[1]*[1])", 3.52, 3.6)
        # g2.SetParLimits(1, 3.4, 3.52)
        # g3.SetParLimits(1, 3.52, 3.7)
        
        g2.SetParLimits(2, 3.485, 3.52)
        g3.SetParLimits(2, 3.53, 3.57)
        g2.SetParameter(0, eephoton1_deltaM.GetMaximum())
        g3.SetParameter(0, eephoton2_deltaM.GetMaximum())
        # g2.SetParLimits(0, 0.85*eephoton1_deltaM.GetMaximum(), 1.2*eephoton1_deltaM.GetMaximum())
        # g3.SetParLimits(0, 0.85*eephoton2_deltaM.GetMaximum(), 1.2*eephoton2_deltaM.GetMaximum())
        g2.SetParameter(2, 3.510)
        g3.SetParameter(2, 3.556)
        g2.SetParameter(1, 0.00084)
        g3.SetParameter(1, 0.03)
        g2.SetParLimits(1, 0.0005, 0.002)
        # g3.SetParLimits(2, 3.53, 3.57)
        total = ROOT.TF1("total", "lorentz1+lorentz2", 3.4, 3.7)
        # total = ROOT.TF1("total", "fGaussExp1+fGaussExp2", 3.44, 3.58)
    elif option == 4: 
        g2 = ROOT.TF1("fGaussExp1",
                 "(x<[1]) * ([0]*(TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2)) + TMath::Exp((x-[1])/[3])*(1. - TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2))))) + \
                 (x>=[1]) * ([0]*TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2)))", 3.44, 3.53)
        g3 = ROOT.TF1("fGaussExp2",
                 "(x<[1]) * ([0]*(TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2)) + TMath::Exp((x-[1])/[3])*(1. - TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2))))) + \
                 (x>=[1]) * ([0]*TMath::Exp(-0.5*TMath::Power((x-[1])/[2],2)))", 3.5, 3.58)
        g2.SetParLimits(1, 3.485, 3.52)
        g3.SetParLimits(1, 3.53, 3.57)
        g2.SetParameter(0, eephoton1_data_deltaM.GetMaximum())
        g3.SetParameter(0, eephoton2_data_deltaM.GetMaximum())
        g2.SetParLimits(0, 0.85*eephoton1_data_deltaM.GetMaximum(), 1.2*eephoton1_data_deltaM.GetMaximum())
        g3.SetParLimits(0, 0.85*eephoton2_data_deltaM.GetMaximum(), 1.2*eephoton2_data_deltaM.GetMaximum())
        g2.SetParameter(1, 3.510)
        g3.SetParameter(1, 3.556)
        g2.SetParameter(2, 0.03)
        g3.SetParameter(2, 0.03)
        g2.SetParameter(3, 0.03)
        g3.SetParameter(3, 0.03)
        total = ROOT.TF1("total", "fGaussExp1+fGaussExp2", 3.44, 3.57)
    elif option == 5: 
        # g2 = ROOT.TF1("lorentz1", "(0.5*[0]*[1]/TMath::Pi()) / TMath::Max(1.e-10,(x-[2])*(x-[2])+ .25*[1]*[1])", 3.4, 3.6)
        # g3 = ROOT.TF1("lorentz2", "(0.5*[0]*[1]/TMath::Pi()) / TMath::Max(1.e-10,(x-[2])*(x-[2])+ .25*[1]*[1])", 3.52, 3.6)
        g2 = ROOT.TF1("lorentz1", "(0.5*[0]*[1]/TMath::Pi()) / ((x-[2])*(x-[2])+ .25*[1]*[1])", 3.485, 3.56)
        g3 = ROOT.TF1("lorentz2", "(0.5*[0]*[1]/TMath::Pi()) / ((x-[2])*(x-[2])+ .25*[1]*[1])", 3.52, 3.57)
        # g2.SetParLimits(1, 3.4, 3.52)
        # g3.SetParLimits(1, 3.52, 3.7)
        
        g2.SetParLimits(2, 3.485, 3.52)
        g3.SetParLimits(2, 3.53, 3.57)
        g2.SetParameter(0, eephoton1_data_deltaM.GetMaximum()/2/3.5/3.1416)
        g3.SetParameter(0, eephoton2_data_deltaM.GetMaximum()/2/3.5/3.1415)
        # peakmin = 0.95*eephoton1_data_deltaM.GetMaximum()/2/3.5/3.1416
        # peakmax = 1.1*eephoton1_data_deltaM.GetMaximum()/2/3.5/3.1416
        # g2.SetParLimits(0, peakmin, peakmax )
        # g3.SetParLimits(0, 0.85*eephoton2_deltaM.GetMaximum(), 1.2*eephoton2_deltaM.GetMaximum())
        g2.SetParameter(2, 3.510)
        g3.SetParameter(2, 3.556)
        # g2.SetParameter(1, 0.03)
        # g3.SetParameter(1, 0.03)
        total = ROOT.TF1("total", "lorentz1+lorentz2", 3.49, 3.585)
    
    elif option == 6: 
        g2 = ROOT.TF1("lorentz1brems",
                 "(x<[2]) * ([0]*((0.5*[1]/TMath::Pi()) / ((x-[2])*(x-[2])+ .25*[1]*[1])) + TMath::Exp((x-[2])/[3])*(1. - ((0.5*[1]/TMath::Pi()) / ((x-[2])*(x-[2])+ .25*[1]*[1])))) + \
                 (x>=[2]) * ([0]*((0.5*[1]/TMath::Pi()) / ((x-[2])*(x-[2])+ .25*[1]*[1])))", 3.44, 3.53)
        g3 = ROOT.TF1("lorentz2brems",
                 "(x<[2]) * ([0]*((0.5*[1]/TMath::Pi()) / ((x-[2])*(x-[2])+ .25*[1]*[1])) + TMath::Exp((x-[2])/[3])*(1. - ((0.5*[1]/TMath::Pi()) / ((x-[2])*(x-[2])+ .25*[1]*[1])))) + \
                 (x>=[2]) * ([0]*((0.5*[1]/TMath::Pi()) / ((x-[2])*(x-[2])+ .25*[1]*[1])))", 3.5, 3.58)
        # g2 = ROOT.TF1("lorentz1", "(0.5*[0]*[1]/TMath::Pi()) / TMath::Max(1.e-10,(x-[2])*(x-[2])+ .25*[1]*[1])", 3.4, 3.6)
        # g3 = ROOT.TF1("lorentz2", "(0.5*[0]*[1]/TMath::Pi()) / TMath::Max(1.e-10,(x-[2])*(x-[2])+ .25*[1]*[1])", 3.52, 3.6)
        # g2 = ROOT.TF1("lorentz1", "(0.5*[0]*[1]/TMath::Pi()) / ((x-[2])*(x-[2])+ .25*[1]*[1])", 3.485, 3.56)
        # g3 = ROOT.TF1("lorentz2", "(0.5*[0]*[1]/TMath::Pi()) / ((x-[2])*(x-[2])+ .25*[1]*[1])", 3.52, 3.57)
        # g2.SetParLimits(1, 3.4, 3.52)
        # g3.SetParLimits(1, 3.52, 3.7)
        
        g2.SetParLimits(2, 3.485, 3.52)
        g3.SetParLimits(2, 3.53, 3.57)
        g2.SetParameter(0, eephoton1_data_deltaM.GetMaximum()/2/3.5/3.1416)
        g3.SetParameter(0, eephoton2_data_deltaM.GetMaximum()/2/3.5/3.1415)
        peakmin = 0.35*eephoton1_data_deltaM.GetMaximum()/2/3.5/3.1416
        peakmax = 1.2*eephoton1_data_deltaM.GetMaximum()/2/3.5/3.1416
        g2.SetParLimits(0, peakmin, peakmax )
        # g3.SetParLimits(0, 0.85*eephoton2_deltaM.GetMaximum(), 1.2*eephoton2_deltaM.GetMaximum())

        g2.SetParameter(2, 3.510)
        g3.SetParameter(2, 3.556)
        g2.SetParameter(3, 0.02)
        g3.SetParameter(3, 0.02)
        g2.SetParLimits(3, 0.04, 0.1)
        g3.SetParLimits(3, 0.04, 0.1)
        # g2.SetParameter(1, 0.03)
        # g3.SetParameter(1, 0.03)
        total = ROOT.TF1("total", "lorentz1brems+lorentz2brems", 3.46, 3.585)


   
        
        
    total.SetLineColor(1)
    g2.SetLineColor(1)
    g3.SetLineColor(1)
    total.SetLineWidth(1)
    g2.SetLineWidth(1)
    g3.SetLineWidth(1)
    
    par = np.zeros(6, dtype=c_double)
    par2 = np.zeros(6, dtype=c_double)
    par3 = np.zeros(8, dtype=c_double)
    par4 = np.zeros(9, dtype=c_double)
    if option == 0:
        eephoton1_deltaM.Fit(g2,"R0")
        eephoton2_deltaM.Fit(g3, "R0+")
        g2.GetParameters(par[0:3])
        g3.GetParameters(par[3:6])
        total.SetParameters(par)
        eephoton12_deltaM.Fit(total, "R+")
        g2.SetParameter(0, total.GetParameter(0))
        g2.SetParameter(1, total.GetParameter(1))
        g2.SetParameter(2, total.GetParameter(2))
        g3.SetParameter(0, total.GetParameter(3))
        g3.SetParameter(1, total.GetParameter(4))
        g3.SetParameter(2, total.GetParameter(5))
        #eephoton12_deltaM.Fit(total, "R+")

        
    if option == 1:
        eephoton12_data_deltaM.Fit(g2,"R")
        eephoton12_data_deltaM.Fit(g3, "R+")
        g2.GetParameters(par[0:3])
        g3.GetParameters(par[3:6])

        total.SetParameters(par)
        eephoton12_data_deltaM.Fit(total, "R+")

        g2.SetParameter(0, total.GetParameter(0))
        g2.SetParameter(1, total.GetParameter(1))
        g2.SetParameter(2, total.GetParameter(2))
        g3.SetParameter(0, total.GetParameter(3))
        g3.SetParameter(1, total.GetParameter(4))
        g3.SetParameter(2, total.GetParameter(5))
        total.GetParameters(par2[0:6])

    if option == 2:
        eephoton1_deltaM.Fit(g2,"R")
        #eephoton12_data_deltaM.Fit(g3, "R+")
        g2.GetParameters(par["Error"])
        #g3.GetParameters(par[3])
        #total.SetParameters(par)
        #eephoton12_data_deltaM.Fit(total, "R+")

        # g2.SetParameter(0, total.GetParameter(0))
        # g2.SetParameter(1, total.GetParameter(1))
        # g2.SetParameter(2, total.GetParameter(2))
        # g3.SetParameter(0, total.GetParameter(3))
        # g3.SetParameter(1, total.GetParameter(4))
        # g3.SetParameter(2, total.GetParameter(5))
    if option == 3: 
        eephoton2_deltaM.Fit(g3,"R")
        g3.GetParameters(par[3])
    if option == 4:
        eephoton1_data_deltaM.Fit(g2,"R0")
        eephoton2_data_deltaM.Fit(g3, "R0+")
        g2.GetParameters(par3[0:4])
        g3.GetParameters(par3[4:8])

        total.SetParameters(par3)
        total.SetParLimits(0, 0.96*eephoton1_data_deltaM.GetMaximum(), 1.2*eephoton1_data_deltaM.GetMaximum())
        total.SetParLimits(1, 3.485, 3.52)
        total.SetParLimits(2, 0.003, 0.007)
        total.SetParLimits(3, 0.004, 0.04)
        total.SetParLimits(4, 0.99*eephoton2_data_deltaM.GetMaximum(), 1.2*eephoton2_data_deltaM.GetMaximum())
        total.SetParLimits(5, 3.53, 3.57)
        total.SetParLimits(6, 0.003, 0.009)
        total.SetParLimits(7, 0.004, 0.04)
        eephoton12_data_deltaM.Fit(total, "R+")

        # g2.SetParameter(0, total.GetParameter(0))
        # g2.SetParameter(1, total.GetParameter(1))
        # g2.SetParameter(2, total.GetParameter(2))
        # g2.SetParameter(3, total.GetParameter(3))
        # g3.SetParameter(0, total.GetParameter(4))
        # g3.SetParameter(1, total.GetParameter(5))
        # g3.SetParameter(2, total.GetParameter(6))
        # g3.SetParameter(3, total.GetParameter(7))
        total.GetParameters(par3[0:8])
        
    if option == 5: 
        eephoton1_data_deltaM.Fit(g2,"R0")
        eephoton2_data_deltaM.Fit(g3, "R0+")
        g2.GetParameters(par[0:3])
        g3.GetParameters(par[3:6])
        total.SetParameters(par)
        # peakmin = 0.98*eephoton12_data_deltaM.GetMaximum()/2/3.2/3.1416
        # peakmax = 1.3*eephoton12_data_deltaM.GetMaximum()/2/3.2/3.1416
        # total.SetParLimits(0, peakmin, peakmax )
        #total.SetParLimits(0, 0.95*eephoton2_data_deltaM.GetMaximum(), 1.2*eephoton2_data_deltaM.GetMaximum())
        total.SetParLimits(5, 3.53, 3.57)
        total.SetParLimits(2, 3.485, 3.52)
        eephoton12_data_deltaM.Fit(total, "R+")
        g2.SetParameter(0, total.GetParameter(0))
        g2.SetParameter(1, total.GetParameter(1))
        g2.SetParameter(2, total.GetParameter(2))
        g3.SetParameter(0, total.GetParameter(3))
        g3.SetParameter(1, total.GetParameter(4))
        g3.SetParameter(2, total.GetParameter(5))
    if option == 6: 
        eephoton1_data_deltaM.Fit(g2,"R0")
        eephoton2_data_deltaM.Fit(g3, "R0+")
        g2.GetParameters(par3[0:4])
        g3.GetParameters(par3[4:8])
        total.SetParameters(par3)
        peakmin = 0.7*eephoton12_data_deltaM.GetMaximum()/2/3.2/3.1416
        peakmax = 1*eephoton12_data_deltaM.GetMaximum()/2/3.2/3.1416
        total.SetParLimits(0, peakmin, peakmax )
        #total.SetParLimits(0, 0.95*eephoton2_data_deltaM.GetMaximum(), 1.2*eephoton2_data_deltaM.GetMaximum())
        total.SetParLimits(6, 3.53, 3.57)
        total.SetParLimits(2, 3.485, 3.52)
        total.SetParLimits(3, 0.02, 0.05)
        total.SetParLimits(7, 0.05, 0.1)
        total.SetParLimits(1, 0.02, 0.024)
        total.SetParLimits(5, 0.016, 0.0185)
        eephoton12_data_deltaM.Fit(total, "R+")

        # g2.SetParameter(0, total.GetParameter(0))
        # g2.SetParameter(1, total.GetParameter(1))
        # g2.SetParameter(2, total.GetParameter(2))
        # g3.SetParameter(0, total.GetParameter(3))
        # g3.SetParameter(1, total.GetParameter(4))
        # g3.SetParameter(2, total.GetParameter(5))
    ymax = 100
  
    if option == 0 : 
        ymax =eephoton12_deltaM.GetMaximum()*2
 
    if option == 1 or option == 4 or option == 5 or option == 6:
        ymax = eephoton12_data_deltaM.GetMaximum()
    if option == 2: 
        ymax = eephoton1_deltaM.GetMaximum()
    if option == 3: 
        ymax = eephoton2_deltaM.GetMaximum()

    
    c1 = TCanvas("DeltaMass","\Delta Mass",0,0,900,900);
    p1 = c1.cd();
    p1.SetPad(0,0.01,0.99,1);
    p1.SetMargin(0.15,0.05,0.12,0.03);
    p1.SetTicks(1,1);
    if option != 1 and option !=4 and option != 5 and option != 6:
        p1.SetLogy() #log scale
    xmin = 3.2
    xmax = 4.2
    if option == 1 or option == 4 or option == 5 or option == 6:
        xmin = 3.4
        xmax = 3.6
    
    if option == 0:
        xmin = 3.4
        xmax = 3.6
    ymin = 0
    if option == 0:
        ymin = 99999


    frame1 = p1.DrawFrame(xmin, ymin,xmax, ymax+(ymax*0.2));
    frame1.GetXaxis().SetTitle("\Deltam + m_{J/\psi}^{PDG} [GeV/c^{2}]");
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
        txt5 = TPaveText(0.6,0.5,0.95,0.57,"NDC");
        txt5.SetFillColor(kWhite);
        txt5.SetFillStyle(0);
        txt5.SetBorderSize(0);
        txt5.SetTextAlign(33);#middle,left
        txt5.SetTextFont(42);#helvetica
        txt5.SetTextSize(0.02);
        txt6 = TPaveText(0.6,0.43,0.95,0.50,"NDC");
        txt6.SetFillColor(kWhite);
        txt6.SetFillStyle(0);
        txt6.SetBorderSize(0);
        txt6.SetTextAlign(33);#middle,left
        txt6.SetTextFont(42);#helvetica
        txt6.SetTextSize(0.02);
     
    if option == 1 or option == 4 or option == 5 or option == 6:
        txt5 = TPaveText(0.15,0.7,0.5,0.77,"NDC");
        txt5.SetFillColor(kWhite);
        txt5.SetFillStyle(0);
        txt5.SetBorderSize(0);
        txt5.SetTextAlign(33);#middle,left
        txt5.SetTextFont(42);#helvetica
        txt5.SetTextSize(0.02);
        txt6 = TPaveText(0.15,0.63,0.5,0.70,"NDC");
        txt6.SetFillColor(kWhite);
        txt6.SetFillStyle(0);
        txt6.SetBorderSize(0);
        txt6.SetTextAlign(33);#middle,left
        txt6.SetTextFont(42);#helvetica
        txt6.SetTextSize(0.02);
        
        

    if option == 0: 
        leg = TLegend(0.17,0.6,.5,0.7);
        eephoton12_deltaM.Draw("Esame")
        #eephoton12_deltaM.Draw("Esame,hist")
        # g2.Draw("same")
        # g3.Draw("same")
        # total.Draw("same")
        eephoton1_deltaM.Draw("Esame")
        eephoton2_deltaM.Draw("Esame")
        leg.AddEntry(eephoton12_deltaM, "\chi_{c1,c2} \\rightarrow \gamma e^{+} e^{-}","LP")
        leg.AddEntry(eephoton1_deltaM, "\chi_{c1} \\rightarrow \gamma e^{+} e^{-}","LP");
        leg.AddEntry(eephoton2_deltaM, "\chi_{c2} \\rightarrow \gamma e^{+} e^{-}","LP");
        # txt5.AddText("\chi_{c1}: M = ("+str(round(par[1],6))+" \pm "+str('{0:.6f}'.format(g2.GetParError(1)))+") GeV/c^{2}  ")
        # #txt5.AddText("\chi_{c1}: M = "+str(round(par[1],3))+" GeV/c^{2}                    ") 
        # txt5.AddText("\chi_{c1}:  \sigma = ("+str(round(par[2]*1000,3))+" \pm "+str('{0:.3f}'.format(g2.GetParError(2)*1000))+") MeV/c^{2}          ")
        
        # #txt5.AddText("\chi_{c1}: \sigma = "+str(round(par[2]*1000,2))+" MeV/c^{2} ")
        # #txt6.AddText("\chi_{c2}: M = "+str(round(par[4],3))+" GeV/c^{2}                    ")
        # txt6.AddText("\chi_{c2}: M = ("+str(round(par[4],6))+" \pm "+str('{0:.6f}'.format(g3.GetParError(1)))+") GeV/c^{2} ")
        # #txt6.AddText("\chi_{c2}: \sigma = "+str(round(par[5]*1000,2))+" MeV/c^{2} ")
        # txt6.AddText("\chi_{c2}:  \sigma = ("+str(round(par[5]*1000,3))+" \pm "+str('{0:.3f}'.format(g3.GetParError(2)*1000))+") MeV/c^{2}          ")


    if option == 1: 
        leg = TLegend(0.17,0.8,.5,0.9);
        # eephoton12_data_deltaM.Draw("Esame")
        eephoton12_data_deltaM.Draw("Esame,hist")
        total.Draw("same")

        eephoton1_data_deltaM.Draw("Esame")
        eephoton2_data_deltaM.Draw("Esame")
        leg.AddEntry(eephoton12_data_deltaM, "\gamma e^{+} e^{-} from \chi_{c1} and \chi_{c2} ","LP")
        leg.AddEntry(eephoton1_data_deltaM, "\gamma e^{+} e^{-} from \chi_{c1} ","LP");
        leg.AddEntry(eephoton2_data_deltaM, "\gamma e^{+} e^{-} from \chi_{c2} ","LP");
        # txt5.AddText("\chi_{c1}: M = ("+str(round(par2[1],4))+" \pm "+str('{0:.4f}'.format(total.GetParError(1)))+") GeV/c^{2} ")
        # txt5.AddText("\chi_{c1}:  \sigma = ("+str(round(par2[2]*1000,1))+" \pm "+str('{0:.1f}'.format(total.GetParError(2)*1000))+") MeV/c^{2}          ")
        # txt6.AddText("\chi_{c2}: M = ("+str(round(par2[4],4))+" \pm "+str('{0:.4f}'.format(total.GetParError(4)))+") GeV/c^{2} ")
        # #txt6.AddText("\chi_{c2}: \sigma = "+str(round(par[5]*1000,2))+" MeV/c^{2} ")
        # txt6.AddText("\chi_{c2}:  \sigma = ("+str(round(par2[5]*1000,1))+" \pm "+str('{0:.1f}'.format(total.GetParError(5)*1000))+") MeV/c^{2}          ")
    
    if option == 4 or option == 5 or option == 6: 
        leg = TLegend(0.17,0.8,.5,0.9);
        # eephoton12_data_deltaM.Draw("Esame")
        eephoton12_data_deltaM.Draw("Esame,hist")
        # total.SetLineColor(2)
        total.Draw("same")
        # g2.Draw("same")
        # g3.Draw("same")
        eephoton1_data_deltaM.Draw("Esame")
        eephoton2_data_deltaM.Draw("Esame")
        leg.AddEntry(eephoton12_data_deltaM, "\chi_{c1,c2} \\rightarrow \gamma e^{+} e^{-}","LP")
        leg.AddEntry(eephoton1_data_deltaM, "\chi_{c1} \\rightarrow \gamma e^{+} e^{-} ","LP");
        leg.AddEntry(eephoton2_data_deltaM, "\chi_{c2} \\rightarrow \gamma e^{+} e^{-}","LP");
        # txt5.AddText("\chi_{c1}: M = ("+str(round(par3[1],4))+" \pm "+str('{0:.4f}'.format(g2.GetParError(1)))+") GeV/c^{2} ")
        # txt5.AddText("\chi_{c1}:  \sigma = ("+str(round(par3[2]*1000,1))+" \pm "+str('{0:.1f}'.format(g2.GetParError(2)*1000))+") MeV/c^{2}           ")
        # txt6.AddText("\chi_{c2}: M = ("+str(round(par3[5],4))+" \pm "+str('{0:.4f}'.format(g3.GetParError(1)))+") GeV/c^{2} ")
        # #txt6.AddText("\chi_{c2}: \sigma = "+str(round(par[5]*1000,2))+" MeV/c^{2} ")
        # txt6.AddText("\chi_{c2}:  \sigma = ("+str(round(par3[6]*1000,1))+" \pm "+str('{0:.1f}'.format(g3.GetParError(2)*1000))+") MeV/c^{2}           ")

    
    if option == 2: 
       leg = TLegend(0.6,0.6,1.,0.7);
       eephoton1_deltaM.Draw("Esame") 
       #g2.Draw("same")
       leg.AddEntry(eephoton1_deltaM, "\gamma e^{+} e^{-} from \chi_{c1} ","LP")
         

    if option == 3: 
       leg = TLegend(0.6,0.6,1.,0.7);
       eephoton2_deltaM.Draw("Esame") 
       #g3.Draw("same")
       leg.AddEntry(eephoton2_deltaM, "\gamma e^{+} e^{-} from \chi_{c2} ","LP")
         
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.03);
    leg.Draw("");
    ROOT.SetOwnership(leg,False);  

    txt5.Draw();
    ROOT.SetOwnership(txt5,False);  
    txt6.Draw();
    ROOT.SetOwnership(txt6,False);  
    if option == 0: 
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
    if option == 0: 
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
    if option == 0: 
        txt3 = TPaveText(0.4,0.75,0.4,0.90,"NDC");
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
    if option == 4: 
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
    if option == 0: 
        txt4 = TPaveText(0.3,0.77,0.38,0.8,"NDC");
        txt4.SetFillColor(kWhite);
        txt4.SetFillStyle(0);
        txt4.SetBorderSize(0);
        txt4.SetTextAlign(33);#middle,left
        txt4.SetTextFont(42);#helvetica
        txt4.SetTextSize(0.03);
        txt4.AddText("MC generated");
        txt4.Draw();
        ROOT.SetOwnership(txt4,False);
    
    if option == 0:
        arrow = TArrow( 3.51069, 200000, 3.51069, 99999, 0.02, '|>' )
        arrow.SetFillStyle( 1001 )
        arrow.Draw()
        arrow2 = TArrow( 3.55617, 200000, 3.55617, 99999, 0.02, '|>' )
        arrow2.SetFillStyle( 1001 )
        arrow2.Draw()
    elif option == 4 or option == 5 or option == 6: 
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
    deltamass_plot("AnalysisResults_chicall_20240224.root", 4, "20240225/plot_fitexpodeltamass_matched.pdf")
    deltamass_plot("AnalysisResults_chicall_20240224.root", 4, "20240225/plot_fitexpodeltamass_matched.svg")
    #deltamass_plot("AnalysisResults_chicall_20240224.root", 6, "20240225/plot_fitdeltamass_matched.pdf")
    # deltamass_plot("AnalysisResults_chicall_20240224.root", 5, "20240225/plot_fitexpodeltamass_matched.svg")
    deltamass_plot("AnalysisResults_chicall_20240224.root", 0, "20240225/plot_fitdeltamass_truth.pdf")
    deltamass_plot("AnalysisResults_chicall_20240224.root", 0, "20240225/plot_fitdeltamass_truth.svg")
# option 0: delta Mass_jpsi truth MC for chic1, chic2 and chic12 together
# option 1: delta Mass_jpsi matched MC for chic1, chic2 and chic12 together 