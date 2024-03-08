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

def photon_plot(filename1, option, plotname):
    rootfile_data = TFile.Open(filename1, "READ");
    list_data = rootfile_data.Get("pcm-qc-mc");
    list_data2 = list_data.Get("Generated");
    list_data3 = list_data.Get("V0")
    list_data4 = list_data3.FindObject("nocut")
    # option 0: Comparing all photons with conversion and reconstructed photons
    # photon efficiency: 
    # option 1: conversion probablilty 
    # option 2: reconstruction probability
    # option 3: total reconstruction probability
    # Conversion points: 
    # option 4: Generated hPhotonRZ
    # option 5: qc RZ photon primary
    # option 6: qc RZ photon primary MC
    if option == 0: 
        photon1 = list_data2.FindObject("hPt_Photon")
        photon1.SetDirectory(0);
        make_common_style(photon1, kFullCircle, 1.0, kGreen+3, 1, 0); #kViolet+2     #kRed      #kCyan+3
        ROOT.SetOwnership(photon1, False);
        photon1.Rebin(5)
        photon1.Scale(1, "width")
        photon2 = list_data2.FindObject("hPt_ConvertedPhoton")
        photon2.SetDirectory(0);
        make_common_style(photon2, kFullCircle, 1.0, kGreen+1, 1, 0);#kMagenta+2    #kOrange+7     #kTeal-1
        ROOT.SetOwnership(photon2, False);
        photon2.Rebin(5)
        photon2.Scale(1, "width")
        photon3 = list_data4.FindObject("hPt_Photon_Primary")
        photon3.SetDirectory(0);
        make_common_style(photon3, kFullCircle, 1.0, kSpring+7, 1, 0);#kMagenta+2    #kOrange+7     #kTeal-1
        ROOT.SetOwnership(photon3, False);
        photon3.Rebin(5)
        photon3.Scale(1, "width")
    if option == 1 or option == 2 or option == 3: 
        photon1 = list_data2.FindObject("hPt_Photon")
        photon1.SetDirectory(0);
        make_common_style(photon1, kFullCircle, 1.0, kGreen+3, 1, 0); #kViolet+2     #kRed      #kCyan+3
        ROOT.SetOwnership(photon1, False);
        photon2 = list_data2.FindObject("hPt_ConvertedPhoton")
        photon2.SetDirectory(0);
        make_common_style(photon2, kFullCircle, 1.0, kGreen+1, 1, 0);#kMagenta+2    #kOrange+7     #kTeal-1
        ROOT.SetOwnership(photon2, False);
        photon3 = list_data4.FindObject("hPt_Photon_Primary")
        photon3.SetDirectory(0);
        make_common_style(photon3, kFullCircle, 1.0, kSpring+7, 1, 0);#kMagenta+2    #kOrange+7     #kTeal-1
        ROOT.SetOwnership(photon3, False);
        print(photon1.GetNbinsX())
        arr_rxy3 = np.array([0, 0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.25,1.5,1.75,2.0,2.5,3,3.5,4.,4.5,5, 6,7, 8])
        h1 = TH1D("pT1", "pT1", 20, arr_rxy3) 
        h2 = TH1D("pT2", "pT2", 20, arr_rxy3) 
        h3 = TH1D("pT3", "pT3", 20, arr_rxy3)
        for ir in range(0, len(arr_rxy3)-1):
            r1 = arr_rxy3[ir];
            r2 = arr_rxy3[ir+1];
            bin_r1 = photon1.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = photon1.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = photon1.IntegralAndError(bin_r1, bin_r2, error, "")
            h1.SetBinContent(ir, content)
            h1.SetBinError(ir, error.value)
            bin_r1_mc = photon2.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_mc = photon2.GetXaxis().FindBin(r2 - 1e-6);
            error_mc = c_double(0.0)
            content_mc = photon2.IntegralAndError(bin_r1_mc, bin_r2_mc, error_mc, "")
            h2.SetBinContent(ir, content_mc)
            h2.SetBinError(ir, error_mc.value)
            bin_r1_3 = photon3.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_3 = photon3.GetXaxis().FindBin(r2 - 1e-6);
            error_3 = c_double(0.0)
            content_3 = photon3.IntegralAndError(bin_r1_3, bin_r2_3, error_3, "")
            h3.SetBinContent(ir, content_3)
            h3.SetBinError(ir, error_3.value)
        
    if option == 13: 
        photon1 = list_data2.FindObject("hPt_Photon")
        photon1.SetDirectory(0);
        make_common_style(photon1, kFullCircle, 1.0, kGreen+3, 1, 0); #kViolet+2     #kRed      #kCyan+3
        ROOT.SetOwnership(photon1, False);
        photon2 = list_data2.FindObject("hPt_ConvertedPhoton")
        photon2.SetDirectory(0);
        make_common_style(photon2, kFullCircle, 1.0, kGreen+1, 1, 0);#kMagenta+2    #kOrange+7     #kTeal-1
        ROOT.SetOwnership(photon2, False);
        photon3 = list_data4.FindObject("hPt_Photon_Primary")
        photon3.SetDirectory(0);
        make_common_style(photon3, kFullCircle, 1.0, kSpring+7, 1, 0);#kMagenta+2    #kOrange+7     #kTeal-1
        ROOT.SetOwnership(photon3, False);
        print(photon1.GetNbinsX())
        arr_rxy3 = np.array([0, 0.1,0.2,0.3,0.4,0.5,0.75,1.0,1.5,2.0,3,4.,4.5])
        arr_rxy4 = np.array([0, 0.1,0.2,0.3,0.4,0.5,0.75,1.0,1.5,2.0,3,4.])
        h1 = TH1D("pT1", "pT1", 11, arr_rxy4) 
        h2 = TH1D("pT2", "pT2", 11, arr_rxy4) 
        h3 = TH1D("pT3", "pT3", 11, arr_rxy4)
        for ir in range(0, len(arr_rxy3)-1):
            r1 = arr_rxy3[ir];
            r2 = arr_rxy3[ir+1];
            bin_r1 = photon1.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2 = photon1.GetXaxis().FindBin(r2 - 1e-6);
            error = c_double(0.0)
            content = photon1.IntegralAndError(bin_r1, bin_r2, error, "")
            h1.SetBinContent(ir, content)
            h1.SetBinError(ir, error.value)
            bin_r1_mc = photon2.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_mc = photon2.GetXaxis().FindBin(r2 - 1e-6);
            error_mc = c_double(0.0)
            content_mc = photon2.IntegralAndError(bin_r1_mc, bin_r2_mc, error_mc, "")
            h2.SetBinContent(ir, content_mc)
            h2.SetBinError(ir, error_mc.value)
            bin_r1_3 = photon3.GetXaxis().FindBin(r1 + 1e-6);
            bin_r2_3 = photon3.GetXaxis().FindBin(r2 - 1e-6);
            error_3 = c_double(0.0)
            content_3 = photon3.IntegralAndError(bin_r1_3, bin_r2_3, error_3, "")
            h3.SetBinContent(ir, content_3)
            h3.SetBinError(ir, error_3.value) 
    
    if option == 4:
        conversionall = list_data2.FindObject("hPhotonRZ")

    if option == 5: 
        conversion = list_data4.FindObject("hRZ_Photon_Primary")

    if option == 6: 
        conversion_mc = list_data4.FindObject("hRZ_Photon_Primary_MC")

    if (option == 0): 
        ymax = photon1.GetMaximum()*1.5
        ymin = 100
        xmin = 0
        xmax = 4.5

        c1 = TCanvas("photon_pT_distribution","Photon pT distribution",0,0,900,900);
        p1 = c1.cd();
        p1.SetPad(0,0.01,1,1);
        p1.SetMargin(0.13,0.05,0.12,0.05);
        p1.SetTicks(1,1);
        p1.SetLogy()
        frame1 = p1.DrawFrame(xmin, ymin, xmax, ymax);
        frame1.GetXaxis().SetTitle("p_{T} [GeV/c]");
        frame1.GetYaxis().SetTitle("counts per GeV/c");
        frame1.GetXaxis().SetTitleSize(0.045);
        frame1.GetYaxis().SetTitleSize(0.045);
        frame1.GetXaxis().SetTitleOffset(1.1);
        frame1.GetYaxis().SetTitleOffset(1.3);
        frame1.GetXaxis().SetLabelSize(0.05);
        frame1.GetYaxis().SetLabelSize(0.05);
        frame1.GetYaxis().SetMaxDigits(3);
        frame1.GetXaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetLabelOffset(0.01);
        ROOT.SetOwnership(frame1,False);

        leg = TLegend(0.6,0.65,0.95,0.75);
        photon1.Draw("Esame")
        photon2.Draw("Esame")
        photon3.Draw("Esame")
        leg.AddEntry(photon1, "all \gamma ","LP");
        leg.AddEntry(photon2, "conversion \gamma ","LP");
        leg.AddEntry(photon3, "reconstructed \gamma", "LP")
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.03);
        leg.Draw("");
        ROOT.SetOwnership(leg,False);  

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

        # txt2 = TPaveText(0.845,0.8,0.83,0.925,"NDC");
        # txt2.SetFillColor(kWhite);
        # txt2.SetFillStyle(0);
        # txt2.SetBorderSize(0);
        # txt2.SetTextAlign(33);#middle,left
        # txt2.SetTextFont(42);#helvetica
        # txt2.SetTextSize(0.03);
        # txt2.AddText("this thesis");
        # txt2.Draw();
        # ROOT.SetOwnership(txt2,False);

        txt3 = TPaveText(0.875,0.8,0.86,0.925,"NDC"); #txt3 = TPaveText(0.9,0.75,0.9,0.90,"NDC");
        txt3.SetFillColor(kWhite);
        txt3.SetFillStyle(0);
        txt3.SetBorderSize(0);
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.03);
        txt3.AddText("pp, #sqrt{s} = 13.6TeV");
        txt3.Draw();
        ROOT.SetOwnership(txt3,False);

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        c1.SaveAs(plotname);
    
    if (option == 1 or option == 2 or option == 3 or option == 13):
        c1 = TCanvas("photon_efficiency","Photon efficiency",0,0,900,900);
        p1 = c1.cd();
        p1.SetPad(0,0.01,1,1);
        p1.SetMargin(0.15,0.05,0.1,0.05);
        p1.SetTicks(1,1);
        #p1.SetLogy()
        leg = TLegend(0.2,0.2,0.5,0.3);
        # rebin_n = 20
        # photon2.Rebin(rebin_n)
        # photon1.Rebin(rebin_n)
        # photon3.Rebin(rebin_n)
        # print(photon1.GetNbinsX())
        # print(photon2.GetNbinsX())
        # print(photon3.GetNbinsX())
        h1eff = h1.Clone("h1effi")
        #h1eff = TH1D("h1effi", "efficiency", 0, 0.0, 20.0)
        h1eff.Sumw2()
        h1eff.SetFillColor(kGreen+1)
        h1eff.SetLineColor(kGreen+1)
        h1eff.SetMarkerColor(kGreen+1)
    
        #Adjust y-axis settings
        y = h1eff.GetYaxis()
        if option == 1: 
            h1eff.Divide(h2, h1, 100, 1, option="B")
            make_common_style(h1eff, kFullCircle, 1.0, kGreen+1, 1, 0);#kMagenta+2    #kOrange+7     #kTeal-1
            ROOT.SetOwnership(h1eff, False);
            y.SetTitle("conversion probability [%]")
            h1eff.SetAxisRange(8,12.5,"y")
        if option == 2: 
            h1eff.Divide(h3, h2, 100, 1, option="B")
            y.SetTitle("reconstruction efficiency \\varepsilon [%]")
        if (option == 3 ): 
            h1eff.Divide(h3, h1, 100, 1, option="B")
            y.SetTitle("\gamma reconstruction efficiency \\varepsilon [%]")
        if (option == 13): 
            h1eff.Divide(h3, h1, 100, 1, option="B")
            y.SetTitle("\gamma reconstruction efficiency \\varepsilon [%]")
        
             
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
        
        h1eff.SetMarkerStyle(kFullCircle)
        if option != 13:
            h1eff.GetXaxis().SetRangeUser(0., 4.5)
            h1eff.Draw("Esame")
        else:
            h1eff.GetXaxis().SetRangeUser(0., 4.5)
            h1eff.Draw("Esame")
        if option == 13:
            txt = TPaveText(0.90,0.55,0.9,0.65,"NDC");
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
        # if option == 13:
        #     txt2 = TPaveText(0.845,0.5,0.83,0.625,"NDC");
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
        if option == 13:
            txt3 = TPaveText(0.845,0.5,0.83,0.625,"NDC");
        else:
            txt3 = TPaveText(0.875,0.8,0.86,0.925,"NDC");
        # if option == 13:
        #     txt3 = TPaveText(0.9,0.45,0.9,0.60,"NDC");
        # else:
        #     txt3 = TPaveText(0.2,0.75,0.4,0.90,"NDC");
        txt3.SetFillColor(kWhite);
        txt3.SetFillStyle(0);
        txt3.SetBorderSize(0);
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.03);
        txt3.AddText("pp, #sqrt{s} = 13.6TeV");
        txt3.Draw();
        ROOT.SetOwnership(txt3,False);

        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        c1.SaveAs(plotname);

    

    if (option == 4 or option == 5 or option == 6):
        c1 = TCanvas("pT_distribution","pT distribution",0,0,900,900);
        p1 = c1.cd();
        p1.SetPad(0.0,0.01,1,1);
        p1.SetMargin(0.13,0.15,0.08,0.05);
        p1.SetTicks(1,1);
        frame1 = p1.DrawFrame(-100, 0, 100, 100);
        frame1.GetXaxis().SetTitle("V_{z} [cm]");
        frame1.GetYaxis().SetTitle("R_{xy} [cm]");
        # frame1.GetZaxis().SetTitle("N_{\gamma}")
        frame1.GetXaxis().SetTitleSize(0.045);
        frame1.GetYaxis().SetTitleSize(0.045);
        # frame1.GetZaxis().SetTitleSize(0.048);
        frame1.GetXaxis().SetTitleOffset(0.9);
        frame1.GetYaxis().SetTitleOffset(1.3);
        # frame1.GetZaxis().SetTitleOffset(0.9);
        frame1.GetXaxis().SetLabelSize(0.05);
        frame1.GetYaxis().SetLabelSize(0.05);
        # frame1.GetZaxis().SetLabelSize(0.05);
        frame1.GetYaxis().SetMaxDigits(3);
        frame1.GetXaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetLabelOffset(0.01);
        frame1.GetZaxis().SetLabelOffset(0.01);
        ROOT.SetOwnership(frame1,False);
        if option == 4:
            conversionall.SetMarkerStyle(20)
            conversionall.GetZaxis().SetTitle("N_{\gamma}")
            conversionall.GetXaxis().SetTitle("V_{z} [cm]")
            conversionall.GetYaxis().SetTitle("R_{xy} [cm]")
            conversionall.GetZaxis().SetTitleOffset(1.35)
            conversionall.SetStats(0)
            conversionall.Draw("COLZ")
        if option == 5:
            conversion.SetMarkerStyle(20)
            conversion.GetZaxis().SetTitle("N_{\gamma}")
            conversion.GetXaxis().SetTitle("V_{z} [cm]")
            conversion.GetYaxis().SetTitle("R_{xy} [cm]")
            conversion.GetZaxis().SetTitleOffset(1.35)
            conversion.SetStats(0)
            conversion.Draw("COLZ")
        if option == 6: 
            conversion_mc.SetMarkerStyle(20)
            conversion_mc.GetZaxis().SetTitle("N_{\gamma}")
            conversion_mc.GetXaxis().SetTitle("V_{z} [cm]")
            conversion_mc.GetYaxis().SetTitle("R_{xy} [cm]")
            conversion_mc.GetZaxis().SetTitleOffset(1.35)
            conversion_mc.SetStats(0)
            conversion_mc.Draw("COLZ")
        if option == 6: 
            txt = TPaveText(0.1,0.22,0.42,0.22,"NDC");
        else: 
            txt = TPaveText(0.1,0.21,0.42,0.21,"NDC");
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
        #     txt2 = TPaveText(0.15,0.21,0.34,0.21,"NDC");
        # txt2.SetFillColor(kWhite);
        # txt2.SetFillStyle(0);
        # txt2.SetBorderSize(0);
        # txt2.SetTextAlign(33);#middle,left
        # txt2.SetTextFont(42);#helvetica
        # txt2.SetTextSize(0.03);
        # txt2.AddText("this thesis");
        # txt2.Draw();
        # ROOT.SetOwnership(txt2,False);

        # if option == 6: 
        #     txt3 = TPaveText(0.09,0.19,0.39,0.19,"NDC");
        # else:
        #     txt3 = TPaveText(0.1,0.17,0.4,0.17,"NDC");
        if option == 6: 
            txt3 = TPaveText(0.09,0.19,0.39,0.19,"NDC");
        else: 
            txt3 = TPaveText(0.1,0.17,0.4,0.17,"NDC");
        txt3.SetFillColor(kWhite);
        txt3.SetFillStyle(0);
        txt3.SetBorderSize(0);
        txt3.SetTextAlign(33);#middle,left
        txt3.SetTextFont(42);#helvetica
        txt3.SetTextSize(0.03);
        txt3.AddText("pp, #sqrt{s} = 13.6TeV");
        txt3.Draw();
        ROOT.SetOwnership(txt3,False);
        if (option == 4):
            ROOT.SetOwnership(txt3,False);
            txt4 = TPaveText(0.25,0.13,0.38,0.13,"NDC");
            txt4.SetFillColor(kWhite);
            txt4.SetFillStyle(0);
            txt4.SetBorderSize(0);
            txt4.SetTextAlign(33);#middle,left
            txt4.SetTextFont(42);#helvetica
            txt4.SetTextSize(0.03);
            txt4.AddText("MC generated");
            txt4.Draw();
            ROOT.SetOwnership(txt4,False);
        elif option == 5:
            txt4 = TPaveText(0.1,0.13,0.4,0.13,"NDC");
            txt4.SetFillColor(kWhite);
            txt4.SetFillStyle(0);
            txt4.SetBorderSize(0);
            txt4.SetTextAlign(33);#middle,left
            txt4.SetTextFont(42);#helvetica
            txt4.SetTextSize(0.03);
            txt4.AddText("MC reconstructed");
            txt4.Draw();
            ROOT.SetOwnership(txt4,False);
        elif option == 6:
            txt4 = TPaveText(0.11,0.16,0.39,0.16,"NDC");
            txt4.SetFillColor(kWhite);
            txt4.SetFillStyle(0);
            txt4.SetBorderSize(0);
            txt4.SetTextAlign(33);#middle,left
            txt4.SetTextFont(42);#helvetica
            txt4.SetTextSize(0.03);
            txt4.AddText("MC reconstructed");
            txt4.Draw();
            ROOT.SetOwnership(txt4,False);
            txt5 = TPaveText(0.1,0.13,0.43,0.13,"NDC");
            txt5.SetFillColor(kWhite);
            txt5.SetFillStyle(0);
            txt5.SetBorderSize(0);
            txt5.SetTextAlign(33);#middle,left
            txt5.SetTextFont(42);#helvetica
            txt5.SetTextSize(0.03);
            txt5.AddText("(MC true coordinates)");
            txt5.Draw();
            ROOT.SetOwnership(txt5,False);
        if option == 4 or option == 5: 
            line = TLine( -15, 1.9, 15, 1.9)
            line.SetLineStyle( 1 )
            line.SetLineWidth(2)
            line.SetLineColor(2)
            line.Draw()
            line1 = TLine( -15, 2.455, 15,2.455)
            line1.SetLineStyle( 1 )
            line1.SetLineWidth(2)
            line1.SetLineColor(2)
            line1.Draw()
            line2 = TLine( -15, 3.235, 15, 3.235 )
            line2.SetLineStyle(1)
            line2.SetLineWidth(2)
            line2.SetLineColor(2)
            line2.Draw()
            line3 = TLine( -15, 3.995, 15, 3.995 )
            line3.SetLineStyle(1)
            line3.SetLineWidth(2)
            line3.SetLineColor(2)
            line3.Draw()
            line4 = TLine( -30, 19.605, 30, 19.605)
            line4.SetLineStyle(1)
            line4.SetLineWidth(2)
            line4.SetLineColor(2)
            line4.Draw()
            line5 = TLine( -35, 24.545, 35, 24.545)
            line5.SetLineStyle(1)
            line5.SetLineWidth(2)
            line5.SetLineColor(2)
            line5.Draw()
            line15 = TLine( -40, 28.5, 40, 28.5)
            line15.SetLineStyle(1)
            line15.SetLineWidth(2)
            line15.SetLineColor(2)
            line15.Draw()
            line6 = TLine( -47, 34.385, 47, 34.385)
            line6.SetLineStyle(1)
            line6.SetLineWidth(2)
            line6.SetLineColor(2)
            line6.Draw()
            line7 = TLine( -52, 39.335, 52, 39.335)
            line7.SetLineStyle(1)
            line7.SetLineWidth(2)
            line7.SetLineColor(2)
            line7.Draw()
            line17 = TLine( -55, 43.3, 55, 43.3)
            line17.SetLineStyle(1)
            line17.SetLineWidth(2)
            line17.SetLineColor(2)
            line17.Draw()
            line8 = TLine( -60, 45.5, 60, 45.5)
            line8.SetLineStyle(1)
            line8.SetLineWidth(2)
            line8.SetLineColor(2)
            line8.Draw()
            line9 = TLine( -65, 50.2, 65, 50.2)
            line9.SetLineStyle(1)
            line9.SetLineWidth(2)
            line9.SetLineColor(2)
            line9.Draw()
            line10 = TLine( -70, 54.5, 70, 54.5)
            line10.SetLineStyle(1)
            line10.SetLineWidth(2)
            line10.SetLineColor(2)
            line10.Draw()
            line11 = TLine( -75, 60.65, 75, 60.65) #???
            line11.SetLineStyle(1)
            line11.SetLineWidth(2)
            line11.SetLineColor(2)
            line11.Draw()
            line12 = TLine( -90, 78.8, 90, 78.8) #??
            line12.SetLineStyle(1)
            line12.SetLineWidth(2)
            line12.SetLineColor(2)
            line12.Draw()
        if option == 6:
            line = TLine( -20, 1.9, 20, 1.9)
            line.SetLineStyle( 1 )
            line.SetLineWidth(2)
            line.SetLineColor(2)
            line.Draw()
            line1 = TLine( -20, 2.455, 20, 2.455)
            line1.SetLineStyle( 1 )
            line1.SetLineWidth(2)
            line1.SetLineColor(2)
            line1.Draw()
            line2 = TLine( -20, 3.235, 20, 3.235 )
            line2.SetLineStyle(1)
            line2.SetLineWidth(2)
            line2.SetLineColor(2)
            line2.Draw()
            line3 = TLine( -20, 3.995, 20, 3.995 )
            line3.SetLineStyle(1)
            line3.SetLineWidth(2)
            line3.SetLineColor(2)
            line3.Draw()
            line4 = TLine( -35, 19.605, 35, 19.605)
            line4.SetLineStyle(1)
            line4.SetLineWidth(2)
            line4.SetLineColor(2)
            line4.Draw()
            line5 = TLine( -45, 24.545, 45, 24.545)
            line5.SetLineStyle(1)
            line5.SetLineWidth(2)
            line5.SetLineColor(2)
            line5.Draw()
            line15 = TLine( -50, 28.5, 50, 28.5)
            line15.SetLineStyle(1)
            line15.SetLineWidth(2)
            line15.SetLineColor(2)
            line15.Draw()
            line6 = TLine( -55, 34.385, 55, 34.385)
            line6.SetLineStyle(1)
            line6.SetLineWidth(2)
            line6.SetLineColor(2)
            line6.Draw()
            line7 = TLine( -60, 39.335, 60, 39.335)
            line7.SetLineStyle(1)
            line7.SetLineWidth(2)
            line7.SetLineColor(2)
            line7.Draw()
            line17 = TLine( -65, 43.3, 65, 43.3)
            line17.SetLineStyle(1)
            line17.SetLineWidth(2)
            line17.SetLineColor(2)
            line17.Draw()
            line8 = TLine( -70, 45.5, 70, 45.5)
            line8.SetLineStyle(1)
            line8.SetLineWidth(2)
            line8.SetLineColor(2)
            line8.Draw()
            line9 = TLine( -75, 50.2, 75, 50.2)
            line9.SetLineStyle(1)
            line9.SetLineWidth(2)
            line9.SetLineColor(2)
            line9.Draw()
            line10 = TLine( -80, 54.5, 80, 54.5)
            line10.SetLineStyle(1)
            line10.SetLineWidth(2)
            line10.SetLineColor(2)
            line10.Draw()
            line11 = TLine( -85, 60.65, 85, 60.65) #???
            line11.SetLineStyle(1)
            line11.SetLineWidth(2)
            line11.SetLineColor(2)
            line11.Draw()
            line12 = TLine( -95, 78.8, 95, 78.8) #??
            line12.SetLineStyle(1)
            line12.SetLineWidth(2)
            line12.SetLineColor(2)
            line12.Draw()



        c1.Modified();
        c1.Update();
        ROOT.SetOwnership(c1,False);
        c1.SaveAs(plotname);


        
            
            
        


if __name__ == "__main__":
    filename = "AnalysisResults_chicall_20240224.root"
    #photon_plot(filename,3, "20240215/plot_totalefficiency_photons.pdf")
    photon_plot(filename,0, "20240305/plot_pT_photons.svg")
    photon_plot(filename,0, "20240305/plot_pT_photons.pdf")
    photon_plot(filename,1, "20240305/plot_photon_conversionprobability.svg")
    photon_plot(filename,1, "20240305/plot_photon_conversionprobability.pdf")
    # photon_plot(filename,2, "20240225/plot_photon_.svg")
    # photon_plot(filename,2, "20240225/plot_.pdf")
    # photon_plot(filename,3, "20240225/plot_photon_reconstructionefficiency.svg")
    # photon_plot(filename,3, "20240225/plot_photon_reconstructionefficiency.pdf")
    photon_plot(filename,4, "20240305/plot_photonconversion_mctruth.png")
    photon_plot(filename,4, "20240305/plot_photonconversion_mctruth.pdf")
    photon_plot(filename,5, "20240305/plot_photonconversion_reconstructed.svg")
    photon_plot(filename,5, "20240305/plot_photonconversion_reconstructed.pdf")
    photon_plot(filename,6, "20240305/plot_photonconversion_mcmatched.svg")
    photon_plot(filename,6, "20240305/plot_photonconversion_mcmatched.pdf")
