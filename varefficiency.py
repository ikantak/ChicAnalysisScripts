import ROOT
import ctypes
from ctypes import *
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency
from ROOT import gStyle, gROOT, gSystem
from ROOT import kWhite, kBlack, kRed, kGreen, kBlue, kYellow, kMagenta, kCyan, kOrange, kAzure, kSpring, kPink, kViolet, kTeal
from ROOT import kFullCircle, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullStar, kFullCross, kFullDiamond, kOpenSquare, kOpenTriangleUp, kOpenCircle, kFullCrossX
import numpy as np
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

def cutefficiency_plot(filename, option, plotname):
    rootfile_data = TFile.Open(filename, "READ");
    list_data = rootfile_data.Get("analysis-dilepton-photon");
    list_data2 = list_data.Get("output");
    if option == 0: 
        data_ee_jpsi = list_data2.FindObject("DileptonsSelected_cut_matchedMC_eeFromJpsi")
        pt_ee_jpsi = data_ee_jpsi.FindObject("Pt")
        
        arr_rxy = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.,2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22])
        h1 = TH1D("pT1", "pT1", 26, arr_rxy) 
        for ir in range(0, len(arr_rxy)-1):
            r1 = arr_rxy[ir];
            r2 = arr_rxy[ir+1];
            bin_r1 = pt_ee_jpsi.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = pt_ee_jpsi.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = pt_ee_jpsi.IntegralAndError(bin_r1, bin_r2, error, "")
            h1.SetBinContent(ir, content)
            h1.SetBinError(ir, error.value)

        mc_ee_jpsi = list_data2.FindObject("MCTruthGenPair_cut_eeFromJpsi")
        pt_mc_ee_jpsi = mc_ee_jpsi.FindObject("Pt")
        h2 = TH1D("pT2", "pT2", 26, arr_rxy) 
        for ir in range(0, len(arr_rxy)-1):
            r1 = arr_rxy[ir];
            r2 = arr_rxy[ir+1];
            bin_r1 = pt_mc_ee_jpsi.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = pt_mc_ee_jpsi.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = pt_mc_ee_jpsi.IntegralAndError(bin_r1, bin_r2, error, "")
            h2.SetBinContent(ir, content)
            h2.SetBinError(ir, error.value)
        
    if option == 1: 
        data_ee_chic1 = list_data2.FindObject("DileptonsSelected_cut_matchedMC_eeFromJpsiFromChic1")
        pt_ee_chic1 = data_ee_chic1.FindObject("Pt")
        pt_ee_chic1.Rebin(20)
        mc_ee_chic1 = list_data2.FindObject("MCTruthGenPair_cut_eeFromJpsiFromChic1")
        pt_mc_ee_chic1 = mc_ee_chic1.FindObject("Pt")
        pt_mc_ee_chic1.Rebin(20)
    if option == 2: 
        data_ee_chic2 = list_data2.FindObject("DileptonsSelected_cut_matchedMC_eeFromJpsiFromChic2")
        pt_ee_chic2 = data_ee_chic2.FindObject("Pt")
        pt_ee_chic2.Rebin(20)
        mc_ee_chic2 = list_data2.FindObject("MCTruthGenPair_cut_eeFromJpsiFromChic2")
        pt_mc_ee_chic2 = mc_ee_chic2.FindObject("Pt")
        pt_mc_ee_chic2.Rebin(20)
    if option ==  11:
        data_ee_chic1 = list_data2.FindObject("DileptonsSelected_cut_matchedMC_eeFromJpsiFromChic1")
        pt_ee_chic1 = data_ee_chic1.FindObject("Pt")
        mc_ee_chic1 = list_data2.FindObject("MCTruthGenPair_cut_eeFromJpsiFromChic1")
        pt_mc_ee_chic1 = mc_ee_chic1.FindObject("Pt")
        arr_rxy = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.,2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22])
        h1 = TH1D("pT1", "pT1", 26, arr_rxy) 
        h2 = TH1D("pT2", "pT2", 26, arr_rxy) 
        for ir in range(0, len(arr_rxy)-1):
            r1 = arr_rxy[ir];
            r2 = arr_rxy[ir+1];
            bin_r1 = pt_ee_chic1.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = pt_ee_chic1.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = pt_ee_chic1.IntegralAndError(bin_r1, bin_r2, error, "")
            h1.SetBinContent(ir, content)
            h1.SetBinError(ir, error.value)
            bin_r1_mc = pt_mc_ee_chic1.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_mc = pt_mc_ee_chic1.GetXaxis().FindBin(r2 - 1e-6);
            error_mc = c_double(0.0)
            content_mc = pt_mc_ee_chic1.IntegralAndError(bin_r1_mc, bin_r2_mc, error_mc, "")
            h2.SetBinContent(ir, content_mc)
            h2.SetBinError(ir, error_mc.value)
        
        
            
        data_ee_chic2 = list_data2.FindObject("DileptonsSelected_cut_matchedMC_eeFromJpsiFromChic2")
        pt_ee_chic2 = data_ee_chic2.FindObject("Pt")
        mc_ee_chic2 = list_data2.FindObject("MCTruthGenPair_cut_eeFromJpsiFromChic2")
        pt_mc_ee_chic2 = mc_ee_chic2.FindObject("Pt")
        h3 = TH1D("pT1", "pT1", 26, arr_rxy) 
        h4 = TH1D("pT2", "pT2", 26, arr_rxy) 
        for ir in range(0, len(arr_rxy)-1):
            r1 = arr_rxy[ir];
            r2 = arr_rxy[ir+1];
            bin_r1 = pt_ee_chic2.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = pt_ee_chic2.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = pt_ee_chic2.IntegralAndError(bin_r1, bin_r2, error, "")
            h3.SetBinContent(ir, content)
            h3.SetBinError(ir, error.value)
            bin_r1_mc = pt_mc_ee_chic2.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_mc = pt_mc_ee_chic2.GetXaxis().FindBin(r2 - 1e-6);
            error_mc = c_double(0.0)
            content_mc = pt_mc_ee_chic2.IntegralAndError(bin_r1_mc, bin_r2_mc, error_mc, "")
            h4.SetBinContent(ir, content_mc)
            h4.SetBinError(ir, error_mc.value)

    if (option == 20): 
        data_photon_chic1 = list_data2.FindObject("Selected_cut_matchedMC_PhotonFromChic1")
        pt_photon_chic1 = data_photon_chic1.FindObject("Pt_Photon")
        mc_photon_chic1 = list_data2.FindObject("MCTruthGen_cut_PhotonFromChic1")
        pt_mc_photon_chic1 = mc_photon_chic1.FindObject("PtMC_photon")
        arr_rxy3 = np.array([0, 0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.25,2.5,2.75,3,3.25,3.5,3.75,4.,4.5,5])
        h1 = TH1D("pT1", "pT1", 22, arr_rxy3) 
        h2 = TH1D("pT2", "pT2", 22, arr_rxy3) 
        for ir in range(0, len(arr_rxy3)-1):
            r1 = arr_rxy3[ir];
            r2 = arr_rxy3[ir+1];
            bin_r1 = pt_photon_chic1.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = pt_photon_chic1.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = pt_photon_chic1.IntegralAndError(bin_r1, bin_r2, error, "")
            h1.SetBinContent(ir, content)
            h1.SetBinError(ir, error.value)
            bin_r1_mc = pt_mc_photon_chic1.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_mc = pt_mc_photon_chic1.GetXaxis().FindBin(r2 - 1e-6);
            error_mc = c_double(0.0)
            content_mc = pt_mc_photon_chic1.IntegralAndError(bin_r1_mc, bin_r2_mc, error_mc, "")
            h2.SetBinContent(ir, content_mc)
            h2.SetBinError(ir, error_mc.value)
        data_photon_chic2 = list_data2.FindObject("Selected_cut_matchedMC_PhotonFromChic2")
        pt_photon_chic2 = data_photon_chic2.FindObject("Pt_Photon")
        mc_photon_chic2 = list_data2.FindObject("MCTruthGen_cut_PhotonFromChic2")
        pt_mc_photon_chic2 = mc_photon_chic2.FindObject("PtMC_photon")
        arr_rxy3 = np.array([0, 0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.25,2.5,2.75,3,3.25,3.5,3.75,4.,4.5,5])
        h3 = TH1D("pT3", "pT3", 22, arr_rxy3) 
        h4 = TH1D("pT4", "pT4", 22, arr_rxy3) 
        for ir in range(0, len(arr_rxy3)-1):
            r1 = arr_rxy3[ir];
            r2 = arr_rxy3[ir+1];
            bin_r1 = pt_photon_chic2.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = pt_photon_chic2.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = pt_photon_chic2.IntegralAndError(bin_r1, bin_r2, error, "")
            h3.SetBinContent(ir, content)
            h3.SetBinError(ir, error.value)
            bin_r1_mc = pt_mc_photon_chic2.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_mc = pt_mc_photon_chic2.GetXaxis().FindBin(r2 - 1e-6);
            error_mc = c_double(0.0)
            content_mc = pt_mc_photon_chic2.IntegralAndError(bin_r1_mc, bin_r2_mc, error_mc, "")
            h4.SetBinContent(ir, content_mc)
            h4.SetBinError(ir, error_mc.value)
        #data_photon_chic12
    if option == 5: 
        data_eephoton_chic1 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic1")
        pt_eephoton_chic1 = data_eephoton_chic1.FindObject("Pt_DileptonPhoton")
        mc_eephoton_chic1 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic1")
        pt_mc_eephoton_chic1 = mc_eephoton_chic1.FindObject("Pt_DileptonPhoton")
    if option == 6: 
        data_eephoton_chic2 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic2")
        pt_eephoton_chic2 = data_eephoton_chic2.FindObject("Pt_DileptonPhoton")
        mc_eephoton_chic2 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic2")
        pt_mc_eephoton_chic2 = mc_eephoton_chic2.FindObject("Pt_DileptonPhoton")
    if option == 7: 
        data_eephoton_chic12 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic12")
        pt_eephoton_chic12 = data_eephoton_chic12.FindObject("Pt_DileptonPhoton")
        mc_eephoton_chic12 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic12")
        pt_mc_eephoton_chic12 = mc_eephoton_chic12.FindObject("Pt_DileptonPhoton")
        arr_rxy2 = np.array([0, 0.5, 1, 1.5, 2., 2.5, 3,  4, 5, 6, 7, 8, 9, 10, 12, 14,16, 18])
        h1 = TH1D("pT1", "pT1", 16, arr_rxy2) 
        h2 = TH1D("pT2", "pT2", 16, arr_rxy2) 
        for ir in range(0, len(arr_rxy2)-1):
            r1 = arr_rxy2[ir];
            r2 = arr_rxy2[ir+1];
            bin_r1 = pt_eephoton_chic12.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = pt_eephoton_chic12.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = pt_eephoton_chic12.IntegralAndError(bin_r1, bin_r2, error, "")
            h1.SetBinContent(ir, content)
            h1.SetBinError(ir, error.value)
            bin_r1_mc = pt_mc_eephoton_chic12.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_mc = pt_mc_eephoton_chic12.GetXaxis().FindBin(r2 - 1e-6);
            error_mc = c_double(0.0)
            content_mc = pt_mc_eephoton_chic12.IntegralAndError(bin_r1_mc, bin_r2_mc, error_mc, "")
            h2.SetBinContent(ir, content_mc)
            h2.SetBinError(ir, error_mc.value)

    if option == 13: 
        data_photon_chic1 = list_data2.FindObject("Selected_cut_matchedMC_PhotonFromChic1")
        pt_photon_chic1 = data_photon_chic1.FindObject("Pt_Photon")
        mc_photon_chic1 = list_data2.FindObject("MCTruthGen_cut_PhotonFromChic1")
        pt_mc_photon_chic1 = mc_photon_chic1.FindObject("PtMC_photon")
        data_photon_chic2 = list_data2.FindObject("Selected_cut_matchedMC_PhotonFromChic2")
        pt_photon_chic2 = data_photon_chic2.FindObject("Pt_Photon")
        mc_photon_chic2 = list_data2.FindObject("MCTruthGen_cut_PhotonFromChic2")
        pt_mc_photon_chic2 = mc_photon_chic2.FindObject("PtMC_photon")
        print(pt_photon_chic1.GetNbinsX())
        print(pt_photon_chic2.GetNbinsX())
        print(pt_mc_photon_chic1.GetNbinsX())
        print(pt_mc_photon_chic2.GetNbinsX())
        pt_photon_chic12= pt_photon_chic1.Clone("photonchic12")
        pt_mc_photon_chic12 = pt_mc_photon_chic1.Clone("mcphotonchic12")
        pt_photon_chic12.Add(pt_photon_chic2)
        pt_mc_photon_chic12.Add(pt_mc_photon_chic2)
        arr_rxy3 = np.array([0, 0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.25,2.5,2.75,3,3.25,3.5,3.75,4.,4.5,5])
        h1 = TH1D("pT1", "pT1", 22, arr_rxy3) 
        h2 = TH1D("pT2", "pT2", 22, arr_rxy3) 
        for ir in range(0, len(arr_rxy3)-1):
            r1 = arr_rxy3[ir];
            r2 = arr_rxy3[ir+1];
            bin_r1 = pt_photon_chic12.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = pt_photon_chic12.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = pt_photon_chic12.IntegralAndError(bin_r1, bin_r2, error, "")
            h1.SetBinContent(ir, content)
            h1.SetBinError(ir, error.value)
            bin_r1_mc = pt_mc_photon_chic12.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_mc = pt_mc_photon_chic12.GetXaxis().FindBin(r2 - 1e-6);
            error_mc = c_double(0.0)
            content_mc = pt_mc_photon_chic12.IntegralAndError(bin_r1_mc, bin_r2_mc, error_mc, "")
            h2.SetBinContent(ir, content_mc)
            h2.SetBinError(ir, error_mc.value)
        # print(pt_photon_chic12.GetNbinsX())
        # print(pt_mc_photon_chic12.GetNbinsX())
        # pt_photon_chic12.Rebin(50)
        # pt_mc_photon_chic12.Rebin(50)
        # print(pt_photon_chic12.GetNbinsX())
        # print(pt_mc_photon_chic12.GetNbinsX())
    if option == 17:
        data_eephoton_chic1 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic1")
        pt_eephoton_chic1 = data_eephoton_chic1.FindObject("Pt_DileptonPhoton")
        mc_eephoton_chic1 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic1")
        pt_mc_eephoton_chic1 = mc_eephoton_chic1.FindObject("Pt_DileptonPhoton")
        data_eephoton_chic2 = list_data2.FindObject("DileptonPhotonInvMass_cut_matchedMC_eePhotonFromChic2")
        pt_eephoton_chic2 = data_eephoton_chic2.FindObject("Pt_DileptonPhoton")
        mc_eephoton_chic2 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic2")
        pt_mc_eephoton_chic2 = mc_eephoton_chic2.FindObject("Pt_DileptonPhoton")
        arr_rxy4 = np.array([0, 0.5, 1, 1.5, 2., 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22])
        h1 = TH1D("pT1", "pT1", 17, arr_rxy4) 
        h2 = TH1D("pT2", "pT2", 17, arr_rxy4) 
        for ir in range(0, len(arr_rxy4)-1):
            r1 = arr_rxy4[ir];
            r2 = arr_rxy4[ir+1];
            bin_r1 = pt_eephoton_chic1.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = pt_eephoton_chic1.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = pt_eephoton_chic1.IntegralAndError(bin_r1, bin_r2, error, "")
            h1.SetBinContent(ir, content)
            h1.SetBinError(ir, error.value)
            bin_r1_mc = pt_mc_eephoton_chic1.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_mc = pt_mc_eephoton_chic1.GetXaxis().FindBin(r2 - 1e-6);
            error_mc = c_double(0.0)
            content_mc = pt_mc_eephoton_chic1.IntegralAndError(bin_r1_mc, bin_r2_mc, error_mc, "")
            h2.SetBinContent(ir, content_mc)
            h2.SetBinError(ir, error_mc.value)
        arr_rxy4 = np.array([0, 0.5, 1, 1.5, 2., 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22])
        h3 = TH1D("pT3", "pT3", 17, arr_rxy4) 
        h4 = TH1D("pT4", "pT4", 17, arr_rxy4) 
        for ir in range(0, len(arr_rxy4)-1):
            r1 = arr_rxy4[ir];
            r2 = arr_rxy4[ir+1];
            bin_r1 = pt_eephoton_chic2.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = pt_eephoton_chic2.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = pt_eephoton_chic2.IntegralAndError(bin_r1, bin_r2, error, "")
            h3.SetBinContent(ir, content)
            h3.SetBinError(ir, error.value)
            bin_r1_mc = pt_mc_eephoton_chic2.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_mc = pt_mc_eephoton_chic2.GetXaxis().FindBin(r2 - 1e-6);
            error_mc = c_double(0.0)
            content_mc = pt_mc_eephoton_chic2.IntegralAndError(bin_r1_mc, bin_r2_mc, error_mc, "")
            h4.SetBinContent(ir, content_mc)
            h4.SetBinError(ir, error_mc.value)
    


    c2 = TCanvas("test","test",0,0,900,900)
    p1 = c2.cd()
    p1.SetPad(0,0.01,1,1)
    p1.SetMargin(0.15,0.05,0.1,0.05)
    p1.SetTicks(1,1)
    
    arr_rxy = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.,2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20])
    arr_rxy2 = np.array([0, 0.5, 1, 1.5, 2., 2.5, 3,  4, 5, 6, 7, 8, 9, 10, 12, 14,16])
    arr_rxy3 = np.array([0, 0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.25,2.5,2.75,3,3.25,3.5,3.75,4.,4.5])
    arr_rxy4 = np.array([0, 0.5, 1, 1.5, 2., 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20])
    if option == 0: 
        h1eff = TH1D("h1effi", "efficiency", 26, arr_rxy)
        h1.Sumw2()
        h2.Sumw2()
        h1eff.Divide(h1, h2, 100, 1, option="B")
    if option == 1: 
        h1eff = TH1D("h1effi", "efficiency", 100, 0, 20)
        pt_ee_chic1.Sumw2()
        pt_mc_ee_chic1.Sumw2()
        h1eff.Divide(pt_ee_chic1, pt_mc_ee_chic1, 100, 1, option="B")
    if option == 2: 
        h1eff = TH1D("h1effi", "efficiency", 100, 0, 20)
        pt_ee_chic2.Sumw2()
        pt_mc_ee_chic2.Sumw2()
        h1eff.Divide(pt_ee_chic2, pt_mc_ee_chic2, 100, 1, option="B")
    if option == 3: 
        h1eff.Divide(pt_photon_chic1, pt_mc_photon_chic1, 100, 1, option="B")
    if option == 4: 
        h1eff.Divide(pt_photon_chic2, pt_mc_photon_chic2, 100, 1, option="B")
    if option == 5: 
        h1eff.Divide(pt_eephoton_chic1, pt_mc_eephoton_chic1, 100, 1, option="B")
    if option == 6: 
        h1eff.Divide(pt_eephoton_chic2, pt_mc_eephoton_chic2, 100, 1, option="B")
    if option == 7: 
        h1eff = TH1D("h1effi", "efficiency", 16,arr_rxy2)
        h1.Sumw2()
        h2.Sumw2()
        h1eff.Divide(h1, h2, 100, 1, option="B")
    if option == 13: 
        h1.Sumw2()
        h2.Sumw2()
        h1eff = TH1D("h1effi", "efficiency", 22, arr_rxy3)
        h1eff.Divide(h1, h2, 100, 1, option="B")
    if option == 11: 
        h1eff = TH1D("h1effi", "efficiency", 26, arr_rxy)
        h1.Sumw2()
        h2.Sumw2()
        h1eff.Divide(h1, h2, 100, 1, option="B")
        h2eff = TH1D("h2effi", "efficiency", 26, arr_rxy)
        h3.Sumw2()
        h4.Sumw2()
        h2eff.Divide(h3, h4, 100, 1, option="B")
    if option == 17: 
        h1eff = TH1D("h1effi", "efficiency", 17, arr_rxy4)
        h1.Sumw2()
        h2.Sumw2()
        h1eff.Divide(h1, h2, 100, 1, option="B")
        h2eff = TH1D("h2effi", "efficiency", 17, arr_rxy4)
        h3.Sumw2()
        h4.Sumw2()
        h2eff.Divide(h3, h4, 100, 1, option="B")
    if option == 20: 
        h1eff = TH1D("h1effi", "efficiency", 22, arr_rxy3)
        h1.Sumw2()
        h2.Sumw2()
        h1eff.Divide(h1, h2, 100, 1, option="B")
        h2eff = TH1D("h2effi", "efficiency", 22, arr_rxy3)
        h3.Sumw2()
        h4.Sumw2()
        h2eff.Divide(h3, h4, 100, 1, option="B")

    #h1eff.Sumw2()
    #h1eff.Scale(100)
    h1eff.SetXTitle("p_{T} [GeV/c]")
    h1eff.SetYTitle("Efficiency")
    

    #Adjust y-axis settings
    y = h1eff.GetYaxis()
    y.SetTitleSize(0.048)
    y.SetTitleFont(42)
    y.SetTitleOffset(1.4)
    y.SetLabelFont(42)
    y.SetLabelSize(0.035)

    #adjust x-axis settings
    x = h1eff.GetXaxis()
    x.SetTitle("p_{T} [GeV/c]")
    x.SetTitleSize(0.048)
    x.SetTitleFont(42)
    x.SetTitleOffset(0.9)
    x.SetLabelFont(42)
    x.SetLabelSize(0.035)
    if option == 0: 
        h1eff.SetFillColor(kCyan+2)
        h1eff.SetMarkerColor(kCyan+2)
        h1eff.SetLineColor(kCyan+2)
    elif option == 1: 
        h1eff.SetFillColor(kCyan+1)
        h1eff.SetMarkerColor(kCyan+1)
        h1eff.SetLineColor(kCyan+1)
    elif option == 2: 
        h1eff.SetFillColor(kCyan-9)
        h1eff.SetMarkerColor(kCyan-9)
        h1eff.SetLineColor(kCyan-9)
    elif option == 7:
        h1eff.SetFillColor(kPink+7)
        h1eff.SetMarkerColor(kPink+7)
        h1eff.SetLineColor(kPink+7)
    elif option == 13:
        h1eff.SetFillColor(kGreen+1)
        h1eff.SetMarkerColor(kGreen+1)
        h1eff.SetLineColor(kGreen+1)
        h1eff.SetAxisRange(0,1.75,"y")
    elif option== 11:
        h1eff.SetFillColor(kCyan+1)
        h1eff.SetMarkerColor(kCyan+1)
        h1eff.SetLineColor(kCyan+1)
        h2eff.SetFillColor(kCyan-9)
        h2eff.SetMarkerColor(kCyan-9)
        h2eff.SetLineColor(kCyan-9)
    elif option == 20:
        h1eff.SetFillColor(kGreen+1)
        h1eff.SetMarkerColor(kGreen+1)
        h1eff.SetLineColor(kGreen+1)
        h1eff.SetAxisRange(0,1.75,"y")
        h2eff.SetFillColor(kSpring+7)
        h2eff.SetMarkerColor(kSpring+7)
        h2eff.SetLineColor(kSpring+7)
        h2eff.SetAxisRange(0,1.75,"y")
    elif option == 17: 
        h1eff.SetFillColor(kPink+7)
        h1eff.SetMarkerColor(kPink+7)
        h1eff.SetLineColor(kPink+7)
        h2eff.SetFillColor(kPink+1)
        h2eff.SetMarkerColor(kPink+1)
        h2eff.SetLineColor(kPink+1)
        h1eff.SetAxisRange(0,0.2,"y")
        h2eff.SetAxisRange(0,0.2,"y")
    else:
        h1eff.SetFillColor(kBlue)
        h1eff.SetMarkerColor(kBlue)
        h1eff.SetLineColor(kBlue)
    
    h1eff.SetMarkerStyle(kFullSquare)
    h1eff.Draw("Esame")
    if option == 11 or option == 17 or option == 20:
        h2eff.SetMarkerStyle(kFullSquare)
        h2eff.Draw("Esame")
    
    #change color 

    if option == 0: 
        leg = TLegend(0.6,0.6,1.0,0.75);
        y.SetTitle("\\varepsilon_{e^{+}e^{-}}^{J/\psi} [%]")
        leg.AddEntry(h1eff, "J/\psi \\rightarrow e^{+} e^{-}", "LP")
    if option == 1: 
        leg = TLegend(0.6,0.6,1.0,0.75);
        y.SetTitle("\\varepsilon_{e^{+}e^{-}}^{\chi_{c1}} [%]")
        leg.AddEntry(h1eff, "e^{+} e^{-} from \chi_{c1}", "LP")
    if option == 2: 
        leg = TLegend(0.6,0.6,1.0,0.75)
        y.SetTitle("\\varepsilon_{e^{+}e^{-}}^{\chic_{c2}} [%]")
        leg.AddEntry(h1eff, "e^{+} e^{-} from \chi_{c2}", "LP")

    if option == 11:
        leg = TLegend(0.6,0.6,1.0,0.7)
        leg.AddEntry(h1eff, "e^{+} e^{-} from \chi_{c1}", "LP")
        leg.AddEntry(h2eff, "e^{+} e^{-} from \chi_{c2}", "LP")
        y.SetTitle("\\varepsilon_{e^{+}e^{-}}^{\chi_{c}} [%]")
    if option == 3: 
        leg = TLegend(0.6,0.6,1.0,0.75)
        leg.AddEntry(h1eff, "\gamma from \chi_{c1}", "LP")
        y.SetTitle("\\varepsilon_{\gamma}^{\chi_{c1}} [%]")
    if option == 4: 
        leg = TLegend(0.6,0.6,1.0,0.75)
        leg.AddEntry(h1eff, "\gamma from \chi_{c2}", "LP")
        y.SetTitle("\\varepsilon_{\gamma}^{\chi_{c2}} [%]")
    if option == 5: 
        leg = TLegend(0.6,0.6,1.0,0.75)
        leg.AddEntry(h1eff, "\gamma e^{+} e^{-} from \chi_{c1}", "LP")
        y.SetTitle("\\varepsilon_{\gamma e^{+}e^{-}}^{\chi_{c1}} [%]")
    if option == 6: 
        leg = TLegend(0.6,0.6,1.0,0.75)
        leg.AddEntry(h1eff, "\gamma e^{+} e^{-} from \chi_{c2}", "LP")
        y.SetTitle("\\varepsilon_{\gamma e^{+}e^{-}}^{\chi_{c2}} [%]")
    if option == 7: 
        leg = TLegend(0.35,0.15,0.55,0.2)
        leg.AddEntry(h1eff, "\gamma e^{+} e^{-} from \chi_{c1} and \chi_{c2}", "LP")
        y.SetTitle("\\varepsilon_{\gamma e^{+}e^{-}}^{\chi_{c}} [%]")
    if option == 8: 
        leg = TLegend(0.6,0.6,1.0,0.75)
        leg.AddEntry(h1eff, "e^{+} e^{-} from PC", "LP")
    if option == 13:
        leg = TLegend(0.3,0.15,0.7,0.25)
        leg.AddEntry(h1eff, "\gamma from \chi_{c1} and \chi_{c2}", "LP")
        y.SetTitle("\\varepsilon_{\gamma}^{\chi_{c}} [%]")
    if option == 17: 
        leg = TLegend(0.5,0.15,0.9,0.2)
        leg.AddEntry(h1eff, "\chi_{c1} \\rightarrow \gamma e^{+} e^{-}", "LP")
        leg.AddEntry(h2eff, "\chi_{c2} \\rightarrow \gamma e^{+} e^{-}", "LP")
        y.SetTitle("\\varepsilon_{\gamma e^{+}e^{-}}^{\chi_{c}} [%]")
    if option == 20:
        leg = TLegend(0.2,0.85,0.55,0.9)
        leg.AddEntry(h1eff, "\gamma from \chi_{c1}", "LP")
        leg.AddEntry(h2eff, "\gamma from \chi_{c2}", "LP")
        y.SetTitle("\\varepsilon_{\gamma}^{\chi_{c}} [%]")
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.03);
    leg.Draw("");
    ROOT.SetOwnership(leg,False);
    if (option == 13 or option ==7 or option == 17): 
        txt = TPaveText(0.2,0.85,0.4,0.95,"NDC");
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
    if (option == 13 or option == 7 or option == 17): 
        txt2 = TPaveText(0.2,0.8,0.35,0.925,"NDC");
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
    if (option == 13 or option == 7 or option == 17): 
        txt3 = TPaveText(0.2,0.75,0.4,0.90,"NDC");
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

    c2.Modified();
    c2.Update();
    ROOT.SetOwnership(c2,False);
    c2.SaveAs(plotname);


if __name__ == "__main__":
    filename = "AnalysisResults_chicall_20240224.root"
    cutefficiency_plot(filename, 20, "20240225/plot_varefficiency_photonchic12.svg")
    cutefficiency_plot(filename, 20, "20240225/plot_varefficiency_photonchic12.pdf")
    # cutefficiency_plot(filename, 17, "20240225/plot_varefficiency_eephotonchic12.svg")
    # cutefficiency_plot(filename, 17, "20240225/plot_varefficiency_eephotonchic12.pdf")
    # cutefficiency_plot(filename, 11, "20240225/plot_varefficiency_eechic12.svg")
    # cutefficiency_plot(filename, 11, "20240225/plot_varefficiency_eechic12.pdf")
    # cutefficiency_plot(filename, 0, "20240225/plot_varefficiency_eejpsi.svg")
    # cutefficiency_plot(filename, 0, "20240225/plot_varefficiency_eejpsi.pdf")
       