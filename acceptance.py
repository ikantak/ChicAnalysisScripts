import re
import numpy as np
import datetime
import math

import ROOT
import ctypes
from ctypes import *
from ROOT import TFile, TDirectory, THashList, TH1F, TH1D, TH2F, TH2, TCanvas, TLegend, TPaveText, TPython, TMath, TF1, TLine, TPython, TEfficiency
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

def cutacceptance_plot(filename, option, plotname):
    #Open the input root file
    rootfile_data = TFile.Open(filename, "READ");
    list_data = rootfile_data.Get("analysis-dilepton-photon");
    list_data2 = list_data.Get("output");
    if option == 0: 
        #import the data from the histogram 
        mc_jpsi = list_data2.FindObject("MCTruthGen_cut_Jpsi")
        pt_mc_jpsi = mc_jpsi.FindObject("PtMC")
        mc_ee_jpsi = list_data2.FindObject("MCTruthGenPair_cut_eeFromJpsi")
        pt_mc_ee_jpsi = mc_ee_jpsi.FindObject("Pt")
        #variable binning size
        arr_rxy = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.,2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22])
        #histograms where the variable binning from the input histograms are saved
        h1 = TH1D("pT1", "pT1", 26, arr_rxy) 
        h2 = TH1D("pT2", "pT2", 26, arr_rxy) 
        for ir in range(0, len(arr_rxy)-1):
            #Identifying the bins for the variable binning area
            r1 = arr_rxy[ir];
            r2 = arr_rxy[ir+1];
            bin_r1 = pt_mc_jpsi.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = pt_mc_jpsi.GetXaxis().FindBin(r2 - 1e-6);
            #before error value
            error = c_double(0.0)
            #integrate over the variable bin and calulate value as well as errror
            content = pt_mc_jpsi.IntegralAndError(bin_r1, bin_r2, error, "")
            #Fill the new histograms with value and error
            h2.SetBinContent(ir, content)
            h2.SetBinError(ir, error.value)
            #same for second histogram
            bin_r1_mc = pt_mc_ee_jpsi.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_mc = pt_mc_ee_jpsi.GetXaxis().FindBin(r2 - 1e-6);
            error_mc = c_double(0.0)
            content_mc = pt_mc_ee_jpsi.IntegralAndError(bin_r1_mc, bin_r2_mc, error_mc, "")
            h1.SetBinContent(ir, content_mc)
            h1.SetBinError(ir, error_mc.value)
    # if (option == 1 or option == 3 or option == 5): 
    #     mc_chic1 = list_data2.FindObject("MCTruthGen_cut_Chic1")
    #     pt_mc_chic1 = mc_chic1.FindObject("Pt")
    #     mc_ee_chic1 = list_data2.FindObject("MCTruthGenPair_cut_eeFromJpsiFromChic1")
    #     pt_mc_ee_chic1 = mc_ee_chic1.FindObject("Pt")
    #     mc_photon_chic1 = list_data2.FindObject("MCTruthGen_cut_PhotonFromChic1")
    #     pt_mc_photon_chic1 = mc_photon_chic1.FindObject("Pt")
    #     mc_eephoton_chic1 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic1")
    #     pt_mc_eephoton_chic1 = mc_eephoton_chic1.FindObject("Pt_DileptonPhoton")
    # if (option == 2 or option == 4 or option == 6): 
    #     mc_chic2 = list_data2.FindObject("MCTruthGen_cut_Chic2")
    #     pt_mc_chic2 = mc_chic2.FindObject("Mass_Pt")
    #     mc_ee_chic2 = list_data2.FindObject("MCTruthGenPair_cut_eeFromJpsiFromChic2")
    #     pt_mc_ee_chic2 = mc_ee_chic2.FindObject("Pt")
    #     mc_photon_chic2 = list_data2.FindObject("MCTruthGen_cut_PhotonFromChic2")
    #     pt_mc_photon_chic2 = mc_photon_chic2.FindObject("Pt")
    #     mc_eephoton_chic2 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic2")
    #     pt_mc_eephoton_chic2 = mc_eephoton_chic2.FindObject("Pt_DileptonPhoton")
    if option == 7: 
        mc_chic1 = list_data2.FindObject("MCTruthGen_cut_Chic1")
        pt_mc_chic1 = mc_chic1.FindObject("PtMC")
        mc_chic2 = list_data2.FindObject("MCTruthGen_cut_Chic2")
        pt_mc_chic2 = mc_chic2.FindObject("PtMC")

        mc_eephoton_chic1 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic1")
        pt_mc_eephoton_chic1 = mc_eephoton_chic1.FindObject("Pt_DileptonPhoton")
        mc_eephoton_chic2 = list_data2.FindObject("MCTruthGenTriple_cut_eePhotonFromChic2")
        pt_mc_eephoton_chic2 = mc_eephoton_chic2.FindObject("Pt_DileptonPhoton")
        arr_rxy = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.,2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22])
        h1 = TH1D("pT1", "pT1", 26, arr_rxy) 
        h2 = TH1D("pT2", "pT2", 26, arr_rxy) 
        for ir in range(0, len(arr_rxy)-1):
            r1 = arr_rxy[ir];
            r2 = arr_rxy[ir+1];
            bin_r1 = pt_mc_chic1.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = pt_mc_chic1.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = pt_mc_chic1.IntegralAndError(bin_r1, bin_r2, error, "")
            h2.SetBinContent(ir, content)
            h2.SetBinError(ir, error.value)
            bin_r1_mc = pt_mc_eephoton_chic1.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_mc = pt_mc_eephoton_chic1.GetXaxis().FindBin(r2 - 1e-6);
            error_mc = c_double(0.0)
            content_mc = pt_mc_eephoton_chic1.IntegralAndError(bin_r1_mc, bin_r2_mc, error_mc, "")
            h1.SetBinContent(ir, content_mc)
            h1.SetBinError(ir, error_mc.value)
        arr_rxy = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.,2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22])
        h3 = TH1D("pT3", "pT3", 26, arr_rxy) 
        h4 = TH1D("pT4", "pT4", 26, arr_rxy) 
        for ir in range(0, len(arr_rxy)-1):
            r1 = arr_rxy[ir];
            r2 = arr_rxy[ir+1];
            bin_r1 = pt_mc_chic2.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = pt_mc_chic2.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = pt_mc_chic2.IntegralAndError(bin_r1, bin_r2, error, "")
            h3.SetBinContent(ir, content)
            h3.SetBinError(ir, error.value)
            bin_r1_mc = pt_mc_eephoton_chic2.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_mc = pt_mc_eephoton_chic2.GetXaxis().FindBin(r2 - 1e-6);
            error_mc = c_double(0.0)
            content_mc = pt_mc_eephoton_chic2.IntegralAndError(bin_r1_mc, bin_r2_mc, error_mc, "")
            h4.SetBinContent(ir, content_mc)
            h4.SetBinError(ir, error_mc.value)

    # if option == 8: 
    #     mc_photon = list_data2.FindObject("MCTruthGen_cut_Photon")
    #     pt_mc_photon = mc_photon.FindObject("Mass_Pt")
    #     mc_dielectronFromPC = list_data2.FindObject("MCTruthGenPair_cut_dielectronFromPC")
    #     pt_mc_dielectronFromPC = mc_dielectronFromPC.FindObject("Pt")
    if option == 12: 
        mc_chic1 = list_data2.FindObject("MCTruthGen_cut_Chic1")
        pt_mc_chic1 = mc_chic1.FindObject("Pt")
        mc_ee_chic1 = list_data2.FindObject("MCTruthGenPair_cut_eeFromJpsiFromChic1")
        pt_mc_ee_chic1 = mc_ee_chic1.FindObject("Pt")
        mc_chic2 = list_data2.FindObject("MCTruthGen_cut_Chic2")
        pt_mc_chic2 = mc_chic2.FindObject("Pt")
        mc_ee_chic2 = list_data2.FindObject("MCTruthGenPair_cut_eeFromJpsiFromChic2")
        pt_mc_ee_chic2 = mc_ee_chic2.FindObject("Pt")
        arr_rxy = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.,2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22])
        h1 = TH1D("pT1", "pT1", 26, arr_rxy) 
        h2 = TH1D("pT2", "pT2", 26, arr_rxy) 
        for ir in range(0, len(arr_rxy)-1):
            r1 = arr_rxy[ir];
            r2 = arr_rxy[ir+1];
            bin_r1 = pt_mc_chic1.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = pt_mc_chic1.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = pt_mc_chic1.IntegralAndError(bin_r1, bin_r2, error, "")
            h2.SetBinContent(ir, content)
            h2.SetBinError(ir, error.value)
            bin_r1_mc = pt_mc_ee_chic1.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_mc = pt_mc_ee_chic1.GetXaxis().FindBin(r2 - 1e-6);
            error_mc = c_double(0.0)
            content_mc = pt_mc_ee_chic1.IntegralAndError(bin_r1_mc, bin_r2_mc, error_mc, "")
            h1.SetBinContent(ir, content_mc)
            h1.SetBinError(ir, error_mc.value)
        arr_rxy = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.,2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22])
        h3 = TH1D("pT3", "pT3", 26, arr_rxy) 
        h4 = TH1D("pT4", "pT4", 26, arr_rxy) 
        for ir in range(0, len(arr_rxy)-1):
            r1 = arr_rxy[ir];
            r2 = arr_rxy[ir+1];
            bin_r1 = pt_mc_chic2.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = pt_mc_chic2.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = pt_mc_chic2.IntegralAndError(bin_r1, bin_r2, error, "")
            h3.SetBinContent(ir, content)
            h3.SetBinError(ir, error.value)
            bin_r1_mc = pt_mc_ee_chic2.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_mc = pt_mc_ee_chic2.GetXaxis().FindBin(r2 - 1e-6);
            error_mc = c_double(0.0)
            content_mc = pt_mc_ee_chic2.IntegralAndError(bin_r1_mc, bin_r2_mc, error_mc, "")
            h4.SetBinContent(ir, content_mc)
            h4.SetBinError(ir, error_mc.value)

    #define window for histogram
    c2 = TCanvas("test","test",0,0,900,900)
    p1 = c2.cd()
    #Define size of the histogram and position of the histogram on the window
    p1.SetPad(0,0.01,1,1)
    p1.SetMargin(0.15,0.05,0.1,0.05)
    p1.SetTicks(1,1)
    arr_rxy = np.array([0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2.,2.25, 2.5, 2.75, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 22])
    h1eff = TH1D("h1effi", "efficiency", 26, arr_rxy )
    h2eff = TH1D("h2effi", "efficiency2", 26, arr_rxy)
    #Adjust y-axis settings
    y = h1eff.GetYaxis()
    if option == 0: 
        #Divide MC J/psi->ee by J/psi multiplied by 100 to get percentage
        #option b is for binomial error
        h1eff.Divide(h1, h2, 100, 1, option="B")
        #set title to y-axis
        y.SetTitle("A_{e^{+}e^{-}}^{J/\psi}  [%]")
        #set legend position
        leg = TLegend(0.6,0.6,1.0,0.75);
        #add entry to legend
        leg.AddEntry(h1eff, "J/\psi \\rightarrow e^{+} e^{-}", "LP")
    # if option == 1: 
    #     h1eff.Divide(pt_mc_ee_chic1, pt_mc_chic1, 100, 1, option="B")
    #     y.SetTitle("A_{e^{+}e^{-}}^{\chi_{c1}}  [%]")
    #     leg = TLegend(0.6,0.6,1.0,0.75);
    #     leg.AddEntry(h1eff, "e^{+} e^{-} from \chi_{c1}", "LP")
    # if option == 2: 
    #     h1eff.Divide(pt_mc_ee_chic2, pt_mc_chic2, 100, 1, option="B")
    #     y.SetTitle("A_{e^{+}e^{-}}^{\chi_{c2}}  [%]")
    #     leg = TLegend(0.6,0.6,1.0,0.75)
    #     leg.AddEntry(h1eff, "e^{+} e^{-} from \chi_{c2}", "LP")
    # if option == 3: 
    #     h1eff.Divide(pt_mc_photon_chic1, pt_mc_chic1, 100, 1, option="B")
    #     y.SetTitle("A_{\gamma}^{\chi_{c1}}  [%]")
    #     leg = TLegend(0.6,0.6,1.0,0.75)
    #     leg.AddEntry(h1eff, "\gamma from \chi_{c1}", "LP")
    # if option == 4: 
    #     h1eff.Divide(pt_mc_photon_chic2, pt_mc_chic2, 100, 1, option="B")
    #     y.SetTitle("A_{\gamma}^{\chi_{c2}}  [%]")
    #     leg = TLegend(0.6,0.6,1.0,0.75)
    #     leg.AddEntry(h1eff, "\gamma from \chi_{c2}", "LP")
    # if option == 5:
    #     h1eff.Divide(pt_mc_eephoton_chic1, pt_mc_chic1, 100, 1, option="B")
    #     y.SetTitle("A_{\gamma e^{+}e^{-}}^{\chi_{c1}}  [%]")
    #     leg = TLegend(0.6,0.6,1.0,0.75)
    #     leg.AddEntry(h1eff, "\gamma e^{+} e^{-} from \chi_{c1}", "LP")
    # if option == 6: 
    #     h1eff.Divide(h1, h2, 100, 1, option="B")
    #     y.SetTitle("A_{\gamma e^{+}e^{-}}^{\chi_{c2}}  [%]")
    #     leg = TLegend(0.6,0.6,1.0,0.75)
    #     leg.AddEntry(h1eff, "\gamma e^{+} e^{-} from \chi_{c2}", "LP")
    if option == 7: 
        h1eff.Divide(h1, h2, 100, 1, option="B")
        h2eff.Divide(h4, h3, 100, 1, option ="B")
        y.SetTitle("A_{\gamma e^{+}e^{-}}^{\chi_{c}} [%]")
        leg = TLegend(0.2,0.65,0.50,0.75)
        leg.AddEntry(h1eff, "\chi_{c1} \\rightarrow \gamma e^{+} e^{-}", "LP")
        leg.AddEntry(h2eff, "\chi_{c2} \\rightarrow \gamma e^{+} e^{-}", "LP")
        h1eff.SetAxisRange(20, 85, "y")
        h2eff.SetAxisRange(20, 85, "y")
    # if option == 8: 
    #     h1eff.Divide(pt_mc_dielectronFromPC, pt_mc_photon, 100, 1, option="B")
    #     y.SetTitle("Acceptance  [%]")
    #     leg = TLegend(0.6,0.6,1.0,0.75)
    #     leg.AddEntry(h1eff, "e^{+} e^{-} from PC", "LP")
    if option == 12: 
        h1eff.Divide(h1, h2, 100, 1, option="B")
        h2eff.Divide(h4, h3, 100, 1, option ="B")
        y.SetTitle("A_{e^{+}e^{-}}^{\chi_{c}}  [%]")
        leg = TLegend(0.6,0.2,0.9,0.3)
        leg.AddEntry(h1eff, "e^{+} e^{-} from \chi_{c1}", "LP")
        leg.AddEntry(h2eff, "e^{+} e^{-} from \chi_{c2}", "LP")
        h1eff.SetAxisRange(40, 70, "y")
        h2eff.SetAxisRange(40,70, "y")
    y.SetTitleSize(0.048)
    y.SetTitleFont(42)
    y.SetTitleOffset(1.4)
    y.SetLabelFont(42)
    y.SetLabelSize(0.035)
    #define color of markers
    if option == 0: 
        h1eff.SetFillColor(kCyan+2)
        h1eff.SetMarkerColor(kCyan+2)
        h1eff.SetLineColor(kCyan+2)
    if option == 12:
        h1eff.SetFillColor(kCyan+1)
        h1eff.SetMarkerColor(kCyan+1)
        h1eff.SetLineColor(kCyan+1)
        h2eff.SetFillColor(kCyan-9)
        h2eff.SetMarkerColor(kCyan-9)
        h2eff.SetLineColor(kCyan-9)
    if option == 7: 
        h1eff.SetFillColor(kPink+7)
        h1eff.SetMarkerColor(kPink+7)
        h1eff.SetLineColor(kPink+7)
        h2eff.SetFillColor(kPink+1)
        h2eff.SetMarkerColor(kPink+1)
        h2eff.SetLineColor(kPink+1)
    #settings to x-axis
    x = h1eff.GetXaxis()
    x.SetTitle("p_{T} [GeV/c]")
    x.SetTitleSize(0.048)
    x.SetTitleFont(42)
    x.SetTitleOffset(0.9)
    x.SetLabelFont(42)
    x.SetLabelSize(0.035)
    #marker style
    h1eff.SetMarkerStyle(kFullCross)
    #draw histogram points
    h1eff.Draw("Esame")
    if option == 7 :
        h2eff.SetMarkerStyle(kFullCrossX)
        h2eff.Draw("Esame")
    if option == 12:
        h2eff.SetMarkerStyle(kFullCross)
        h2eff.Draw("Esame")
   
    #settings legend
    leg.SetBorderSize(0);
    leg.SetFillColor(kWhite);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.03);
    leg.Draw("");
    ROOT.SetOwnership(leg,False);
    #Additional text to the histogram
    txt = TPaveText(0.451,0.85,0.451,0.95,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(33);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.03);
    txt.AddText("Simulation this thesis");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    txt3 = TPaveText(0.345,0.8,0.43,0.925,"NDC"); # txt3 = TPaveText(0.4,0.75,0.4,0.90,"NDC");
    txt3.SetFillColor(kWhite);
    txt3.SetFillStyle(0);
    txt3.SetBorderSize(0);
    txt3.SetTextAlign(33);#middle,left
    txt3.SetTextFont(42);#helvetica
    txt3.SetTextSize(0.03);
    txt3.AddText("pp, #sqrt{s} = 13.6TeV");
    txt3.Draw();
    ROOT.SetOwnership(txt3,False);
    #save histogram
    c2.Modified();
    c2.Update();
    ROOT.SetOwnership(c2,False);
    c2.SaveAs(plotname);


if __name__ == "__main__":
    filename = "AnalysisResults_chicall_20240224.root"
    cutacceptance_plot(filename, 0, "20240305/plot_acceptance_eejpsi.pdf")
    cutacceptance_plot(filename, 0, "20240305/plot_acceptance_eejpsi.svg")
    cutacceptance_plot(filename, 7, "20240305/plot_acceptance_eephotonchic.pdf")
    cutacceptance_plot(filename, 7, "20240305/plot_acceptance_eephotonchic.svg")
    cutacceptance_plot(filename, 12, "20240305/plot_acceptance_eechic12.pdf")
    cutacceptance_plot(filename, 12, "20240305/plot_acceptance_eechic12.svg")
    
