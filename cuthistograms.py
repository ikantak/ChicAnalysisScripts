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

def cut_plot(filename1, option, plotname):
    rootfile_data = TFile.Open(filename1, "READ");
    list_data = rootfile_data.Get("pcm-qc-mc");
    list_data2 = list_data.Get("V0Leg");
    list_data3 = list_data.Get("V0")
    list_data4 = list_data3.FindObject("nocut")
    list_data5 = list_data2.FindObject("nocut")
    e_data = rootfile_data.Get("analysis-track-selection")
    e_data_2 = e_data.Get("output")
    e_folder = e_data_2.FindObject("TrackBarrel_jpsiO2MCdebugCuts2")

    cosPA = list_data4.FindObject("hCosPA")
    make_common_style(cosPA, kFullCircle, 1.0, kGreen+1, 1, 0);#kMagenta+2    #kOrange+7     #kTeal-1
    ROOT.SetOwnership(cosPA, False);
    dEdxPCM = list_data5.FindObject("hTPCdEdx")
    hAP = list_data4.FindObject("hAPplot")
    dEdx_e = e_folder.FindObject("TPCnSigEle_pIN")
 
    c1 = TCanvas("pT_distribution","pT distribution",0,0,1500,900);
    p1 = c1.cd();
    p1.SetPad(0.0,0.01,1,1);
    if option == 1: 
        p1.SetMargin(0.13,0.05,0.13,0.05);
    else: 
        p1.SetMargin(0.1,0.12,0.1,0.05);
    p1.SetTicks(1,1);
    xmin = 0
    xmax = 10
    ymin = 0
    ymax = 200
    if option == 1: 
        xmin = 0.985
        ymin = 1
        xmax = 1
        ymax = 1000000
    if option == 2:
        ymin = 40
        ymax = 140
    if option == 3: 
        xmin = -1
        xmax = 1
        ymax = 0.25

    frame1 = p1.DrawFrame(xmin, ymin, xmax, ymax);
    if option == 1: 
        p1.SetLogy()
        frame1.GetXaxis().SetTitle("cos(\\theta_{PA})");
        frame1.GetYaxis().SetTitle("Counts");
        frame1.GetXaxis().SetTitleSize(0.045);
        frame1.GetYaxis().SetTitleSize(0.045);
        frame1.GetXaxis().SetTitleOffset(1.3);
        frame1.GetYaxis().SetTitleOffset(1.4);
        # frame1.GetZaxis().SetTitleOffset(0.9);
        frame1.GetXaxis().SetLabelSize(0.04);
        frame1.GetYaxis().SetLabelSize(0.04);
        # frame1.GetZaxis().SetLabelSize(0.05);
        frame1.GetYaxis().SetMaxDigits(3);
        frame1.GetXaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetLabelOffset(0.01);
        frame1.GetZaxis().SetLabelOffset(0.01);
        ROOT.SetOwnership(frame1,False);
        cosPA.Draw("Esame,hist")
    if option == 2: 
        # frame1.GetXaxis().SetTitle("p_{in} [GeV/c]");
        # frame1.GetYaxis().SetTitle("TPC \\frac{dE}{dx}");
        y = dEdxPCM.GetYaxis()
        y.SetTitle("TPC \\frac{dE}{dx}")
        y.SetTitleSize(0.04)
        y.SetTitleFont(42)
        y.SetTitleOffset(1.)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)
        x = dEdxPCM.GetXaxis()
        x.SetTitle("p_{in} [GeV/c]")
        x.SetTitleSize(0.04)
        x.SetTitleFont(42)
        x.SetTitleOffset(1)
        x.SetLabelFont(42)
        x.SetLabelSize(0.035)
        dEdxPCM.Draw("COLZ")
    if option == 3: 
        # frame1.GetXaxis().SetTitle("\\alpha");
        # frame1.GetYaxis().SetTitle("q_{T} [GeV/c]");
        y = hAP.GetYaxis()
        y.SetTitle("q_{T} [GeV/c]")
        y.SetTitleSize(0.04)
        y.SetTitleFont(42)
        y.SetTitleOffset(1.)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)
        x = hAP.GetXaxis()
        x.SetTitle("\\alpha")
        x.SetTitleSize(0.04)
        x.SetTitleFont(42)
        x.SetTitleOffset(1)
        x.SetLabelFont(42)
        x.SetLabelSize(0.035)
        hAP.Draw("COLZ")
    if option == 4: 
        # frame1.GetXaxis().SetTitle("p_{in} [GeV/c]");
        # frame1.GetYaxis().SetTitle("n \sigma_{e}^{TPC}");
        y = dEdx_e.GetYaxis()
        y.SetTitle("n \sigma_{e}^{TPC}")
        y.SetTitleSize(0.04)
        y.SetTitleFont(42)
        y.SetTitleOffset(1)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)
        x = dEdx_e.GetXaxis()
        x.SetTitle("p_{in} [GeV/c]")
        x.SetTitleSize(0.04)
        x.SetTitleFont(42)
        x.SetTitleOffset(1)
        x.SetLabelFont(42)
        x.SetLabelSize(0.035)
        dEdx_e.Draw("COLZ")
    
    if option == 6: 
        txt = TPaveText(0.08,0.25,0.38,0.25,"NDC");
    else: 
        txt = TPaveText(0.851,0.85,0.851,0.95,"NDC");
    txt.SetFillColor(kWhite);
    txt.SetFillStyle(0);
    txt.SetBorderSize(0);
    txt.SetTextAlign(33);#middle,left
    txt.SetTextFont(42);#helvetica
    txt.SetTextSize(0.03);
    txt.AddText("Simulation this thesis");
    txt.Draw();
    ROOT.SetOwnership(txt,False);

    # if option == 6: 
    #     txt2 = TPaveText(0.14,0.22,0.33,0.22,"NDC");
    # else: 
    #     txt2 = TPaveText(0.81,0.8,0.81, 0.925,"NDC");
    # txt2.SetFillColor(kWhite);
    # txt2.SetFillStyle(0);
    # txt2.SetBorderSize(0);
    # txt2.SetTextAlign(33);#middle,left
    # txt2.SetTextFont(42);#helvetica
    # txt2.SetTextSize(0.03);
    # txt2.AddText("this thesis");
    # txt2.Draw();
    # ROOT.SetOwnership(txt2,False);

    if option == 6: 
        txt3 = TPaveText(0.09,0.19,0.39,0.19,"NDC");
    else:
        txt3 = TPaveText(0.84,0.8,0.84,0.925,"NDC");
    txt3.SetFillColor(kWhite);
    txt3.SetFillStyle(0);
    txt3.SetBorderSize(0);
    txt3.SetTextAlign(33);#middle,left
    txt3.SetTextFont(42);#helvetica
    txt3.SetTextSize(0.03);
    txt3.AddText("pp, #sqrt{s} = 13.6TeV");
    txt3.Draw();
    ROOT.SetOwnership(txt3,False);
    # if (option == 4):
    #     ROOT.SetOwnership(txt3,False);
    #     txt4 = TPaveText(0.25,0.13,0.38,0.13,"NDC");
    #     txt4.SetFillColor(kWhite);
    #     txt4.SetFillStyle(0);
    #     txt4.SetBorderSize(0);
    #     txt4.SetTextAlign(33);#middle,left
    #     txt4.SetTextFont(42);#helvetica
    #     txt4.SetTextSize(0.03);
    #     txt4.AddText("MC generated");
    #     txt4.Draw();
    #     ROOT.SetOwnership(txt4,False);
    # elif option == 5:
    #     txt4 = TPaveText(0.1,0.13,0.4,0.13,"NDC");
    #     txt4.SetFillColor(kWhite);
    #     txt4.SetFillStyle(0);
    #     txt4.SetBorderSize(0);
    #     txt4.SetTextAlign(33);#middle,left
    #     txt4.SetTextFont(42);#helvetica
    #     txt4.SetTextSize(0.03);
    #     txt4.AddText("MC reconstructed");
    #     txt4.Draw();
    #     ROOT.SetOwnership(txt4,False);
    c1.Modified();
    c1.Update();
    ROOT.SetOwnership(c1,False);
    c1.SaveAs(plotname);

        
            
            
        


if __name__ == "__main__":
    filename = "AnalysisResults_chicall_20240224.root"
    cut_plot(filename, 1, "20240305/plot_cosPA.pdf")
    cut_plot(filename, 2, "20240305/plot_dEdx_PCM.pdf")
    cut_plot(filename, 3, "20240305/plot_AP.pdf")
    cut_plot(filename, 4, "20240305/plot_dEdx_ep.pdf")
