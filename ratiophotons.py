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

def ratiophoton_plot(filename1, option, plotname):
    rootfile_data = TFile.Open(filename1, "READ");
    list_data = rootfile_data.Get("pcm-qc-mc");
    list_data2 = list_data.Get("Generated");
    list_data3 = list_data.Get("V0")
    list_data4 = list_data3.FindObject("nocut")
    

    conversionall = list_data2.FindObject("hPhotonRZ")
    conversion = list_data4.FindObject("hRZ_Photon_Primary")
    conversion_mc = list_data4.FindObject("hRZ_Photon_Primary_MC")
    
    rxy_truth = TH1D("R_xy_truth", "R_xy_truth", 100, 0, 100)
    rxy_recon_mc = TH1D("R_xy_reconmc", "R_xy_reconmc", 100, 0, 100 ) 
    rxy_recon = TH1D("R_xy_recon", "R_xy_recon", 100, 0, 100 ) 
    rxy_truth = conversionall.ProjectionY(name="_py", firstxbin=0, lastxbin=-1, option="")
    rxy_recon = conversion.ProjectionY(name="_py", firstxbin=0, lastxbin=-1, option="")
    rxy_recon_mc = conversion_mc.ProjectionY(name="_py", firstxbin=0, lastxbin=-1, option="")
    rxy_truth.Rebin(10)
    # rxy_recon.Rebin(2)
    # rxy_recon_mc.Rebin(2)
    print(rxy_truth.GetNbinsX())
    print(rxy_recon.GetNbinsX())
    print(rxy_recon_mc.GetNbinsX())
        # area = range(1, 101)
    # for value in area:
        # r_truth = conversionall.ProjectionY(name="_px", firstybin=value, lastybin=value+1, option="")
        # r_recon = conversion.ProjectionX(name="_px", firstybin=value, lastybin=value+1, option="")
        # r_recon_mc = conversion_mc.ProjectionX(name="_px", firstybin=value, lastybin=value+1, option="")
        # bin1 = r_truth.GetXaxis().FindBin(0.0000001)
        # bin2 = r_truth.GetXaxis().FindBin(100-0.0000001)
        # error_projection_truth = c_double(0.0)
        # error_projection_recon = c_double(0.0)
        # error_projection_recon_mc = c_double(0.0)
        # integral_truth = r_truth.IntegralAndError(bin1, bin2, error_projection_truth, "")
        # integral_recon = r_recon.IntegralAndError(bin1, bin2, error_projection_recon, "")
        # integral_recon_mc = r_recon_mc.IntegralAndError(bin1, bin2, error_projection_recon_mc, "")
        # rxy_truth.SetBinContent(value, integral_truth)
        # rxy_truth.SetBinError(value, error_projection_truth.value)
        # rxy_recon.SetBinContent(value, integral_recon)
        # rxy_recon.SetBinError(value, error_projection_recon.value)
        # rxy_recon_mc.SetBinContent(value, integral_recon_mc)
        # rxy_recon_mc.SetBinError(value, error_projection_recon_mc.value)

    c1 = TCanvas("pT_distribution","pT distribution",0,0,1500,900);
    p1 = c1.cd();
    p1.SetPad(0.0,0.01,1,1);
    p1.SetMargin(0.1,0.03,0.1,0.03);
    p1.SetTicks(1,1);
    if option == 1: 
        p1.SetLogy()
        make_common_style(rxy_truth, kFullCircle, 1.0, kGreen+3, 1, 0);#kMagenta+2    #kOrange+7     #kTeal-1
        ROOT.SetOwnership(rxy_truth, False);
        make_common_style(rxy_recon, kFullCircle, 1.0, kSpring+7, 1, 0);#kMagenta+2    #kOrange+7     #kTeal-1
        ROOT.SetOwnership(rxy_recon, False);
        make_common_style(rxy_recon_mc, kFullCircle, 1.0, kGreen+1, 1, 0);#kMagenta+2    #kOrange+7     #kTeal-1
        ROOT.SetOwnership(rxy_recon_mc, False);
        # rxy_truth.Scale(1, "width")
        # rxy_recon.Scale(1, "width")
        # rxy_recon_mc.Scale(1, "width")
        xmin = 0
        ymin = 1#min(rxy_truth.GetMinimum(), rxy_recon.GetMinimum(), rxy_recon_mc.GetMinimum())
        ymax = max( rxy_truth.GetMaximum(), rxy_recon.GetMaximum(), rxy_recon_mc.GetMaximum())*3
        xmax = 100
        frame1 = p1.DrawFrame(xmin, ymin, xmax, ymax);
        frame1.GetXaxis().SetTitle("R_{xy} [cm]");
        frame1.GetYaxis().SetTitle("N_{\gamma}");
        frame1.GetXaxis().SetTitleSize(0.045);
        frame1.GetYaxis().SetTitleSize(0.045);
        frame1.GetXaxis().SetTitleOffset(1.1);
        frame1.GetYaxis().SetTitleOffset(1.);
        frame1.GetXaxis().SetLabelSize(0.05);
        frame1.GetYaxis().SetLabelSize(0.05);
        frame1.GetYaxis().SetMaxDigits(3);
        frame1.GetXaxis().SetLabelOffset(0.01);
        frame1.GetYaxis().SetLabelOffset(0.01);
        ROOT.SetOwnership(frame1,False);

        
        rxy_truth.Draw("Esame,hist")
        rxy_recon_mc.Draw("Esame,hist")
        rxy_recon.Draw("Esame,hist")
        line = TLine(1.9, ymin, 1.9, ymax)
        line.SetLineStyle( 2 )
        line.SetLineWidth(2)
        line.SetLineColor(17)
        line.Draw()
        line1 = TLine(2.455, ymin, 2.455, ymax)
        line1.SetLineStyle( 2 )
        line1.SetLineWidth(2)
        line1.SetLineColor(17)
        line1.Draw()
        line2 = TLine(3.235, ymin, 3.235, ymax )
        line2.SetLineStyle(2)
        line2.SetLineWidth(2)
        line2.SetLineColor(17)
        line2.Draw()
        line3 = TLine(3.995, ymin, 3.995, ymax )
        line3.SetLineStyle(2)
        line3.SetLineWidth(2)
        line3.SetLineColor(17)
        line3.Draw()
        line4 = TLine(19.605, ymin, 19.605, ymax)
        line4.SetLineStyle(2)
        line4.SetLineWidth(2)
        line4.SetLineColor(17)
        line4.Draw()
        line5 = TLine(24.545, ymin, 24.545, ymax)
        line5.SetLineStyle(2)
        line5.SetLineWidth(2)
        line5.SetLineColor(17)
        line5.Draw()
        line15 = TLine(28.5, ymin, 28.5, ymax)
        line15.SetLineStyle(2)
        line15.SetLineWidth(2)
        line15.SetLineColor(17)
        line15.Draw()
        line6 = TLine(34.385, ymin, 34.385, ymax)
        line6.SetLineStyle(2)
        line6.SetLineWidth(2)
        line6.SetLineColor(17)
        line6.Draw()
        line7 = TLine(39.335, ymin, 39.335, ymax)
        line7.SetLineStyle(2)
        line7.SetLineWidth(2)
        line7.SetLineColor(17)
        line7.Draw()
        line17 = TLine(43.3, ymin, 43.3, ymax)
        line17.SetLineStyle(2)
        line17.SetLineWidth(2)
        line17.SetLineColor(17)
        line17.Draw()
        line8 = TLine(45.5, ymin, 45.5, ymax)
        line8.SetLineStyle(2)
        line8.SetLineWidth(2)
        line8.SetLineColor(17)
        line8.Draw()
        line9 = TLine(50.2, ymin, 50.2, ymax)
        line9.SetLineStyle(2)
        line9.SetLineWidth(2)
        line9.SetLineColor(17)
        line9.Draw()
        line10 = TLine(54.5, ymin, 54.5, ymax)
        line10.SetLineStyle(2)
        line10.SetLineWidth(2)
        line10.SetLineColor(17)
        line10.Draw()
        line11 = TLine(60.65, ymin, 60.65,ymax) #???
        line11.SetLineStyle(2)
        line11.SetLineWidth(2)
        line11.SetLineColor(17)
        line11.Draw()
        line12 = TLine(78.8, ymin, 78.8, ymax) #??
        line12.SetLineStyle(2)
        line12.SetLineWidth(2)
        line12.SetLineColor(17)
        line12.Draw()
        leg = TLegend(0.4,0.12,0.7,0.22);
        leg.AddEntry(rxy_truth, "MC generated true","LP");
        leg.AddEntry(rxy_recon_mc, "MC reconstructed matched to true coordinates", "LP")
        leg.AddEntry(rxy_recon, "MC reconstructed","LP");
       
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.03);
        leg.Draw("");
        ROOT.SetOwnership(leg,False);  

        txt = TPaveText(0.942,0.85,0.942,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(33);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.03);
        txt.AddText("Simulation this thesis");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        # txt2 = TPaveText(0.89,0.8,0.89, 0.925,"NDC");
        # txt2.SetFillColor(kWhite);
        # txt2.SetFillStyle(0);
        # txt2.SetBorderSize(0);
        # txt2.SetTextAlign(33);#middle,left
        # txt2.SetTextFont(42);#helvetica
        # txt2.SetTextSize(0.03);
        # txt2.AddText("this thesis");
        # txt2.Draw();
        # ROOT.SetOwnership(txt2,False);

        txt3 = TPaveText(0.93,0.8,0.93,0.925,"NDC");#txt3 = TPaveText(0.92,0.75,0.92,0.90,"NDC");
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
    

    if option == 2:     
        rp1 = rxy_truth.Clone("ratio_truth_recon")
        rp2 = rxy_truth.Clone("ratio_truth_recon_mc")
        rp1.Sumw2()
        # rp1.SetFillColor(kGreen+1)
        # rp1.SetLineColor(kGreen+1)
        # rp1.SetMarkerColor(kGreen+1)
        rp2.Sumw2()
        # rp2.SetFillColor(kSpring+7)
        # rp2.SetLineColor(kSpring+7)
        # rp2.SetMarkerColor(kSpring+7)
        rp1.Divide(rxy_recon, rxy_truth, 1, 1, option = "B")
        rp2.Divide(rxy_recon_mc, rxy_truth, 100, 1, option = "B")
        p1.SetLogy()
        y = rp2.GetYaxis()
        y.SetTitle("\\varepsilon_{\gamma} [%]")
        y.SetTitleSize(0.048)
        y.SetTitleFont(42)
        y.SetTitleOffset(1)
        y.SetLabelFont(42)
        y.SetLabelSize(0.035)
        x = rp2.GetXaxis()
        x.SetTitle("R_{xy} [cm]")
        x.SetTitleSize(0.048)
        x.SetTitleFont(42)
        x.SetTitleOffset(0.9)
        x.SetLabelFont(42)
        x.SetLabelSize(0.035)

        # make_common_style(rp1, kFullCircle, 1.0, kSpring+7, 1, 0);#kMagenta+2    #kOrange+7     #kTeal-1
        # ROOT.SetOwnership(rp1, False);
        make_common_style(rp2, kFullCircle, 1.0, kGreen+1, 1, 0);#kMagenta+2    #kOrange+7     #kTeal-1
        ROOT.SetOwnership(rp2, False);
        ymin = min(rp1.GetMinimum(), rp2.GetMinimum())
        ymax = 100
        #rp1.Draw("Esame,hist")
        #
        rp2.GetYaxis().SetRangeUser(0.001, ymax)
        rp2.Draw("Esame,hist")
        
        
        line = TLine(1.9, ymin, 1.9, ymax)
        line.SetLineStyle( 2 )
        line.SetLineWidth(2)
        line.SetLineColor(17)
        line.Draw()
        line1 = TLine(2.455, ymin, 2.455, ymax)
        line1.SetLineStyle( 2 )
        line1.SetLineWidth(2)
        line1.SetLineColor(17)
        line1.Draw()
        line2 = TLine(3.235, ymin, 3.235, ymax )
        line2.SetLineStyle(2)
        line2.SetLineWidth(2)
        line2.SetLineColor(17)
        line2.Draw()
        line3 = TLine(3.995, ymin, 3.995, ymax )
        line3.SetLineStyle(2)
        line3.SetLineWidth(2)
        line3.SetLineColor(17)
        line3.Draw()
        line4 = TLine(19.605, ymin, 19.605, ymax)
        line4.SetLineStyle(2)
        line4.SetLineWidth(2)
        line4.SetLineColor(17)
        line4.Draw()
        line5 = TLine(24.545, ymin, 24.545, ymax)
        line5.SetLineStyle(2)
        line5.SetLineWidth(2)
        line5.SetLineColor(17)
        line5.Draw()
        line15 = TLine(28.5, ymin, 28.5, ymax)
        line15.SetLineStyle(2)
        line15.SetLineWidth(2)
        line15.SetLineColor(17)
        line15.Draw()
        line6 = TLine(34.385, ymin, 34.385, ymax)
        line6.SetLineStyle(2)
        line6.SetLineWidth(2)
        line6.SetLineColor(17)
        line6.Draw()
        line7 = TLine(39.335, ymin, 39.335, ymax)
        line7.SetLineStyle(2)
        line7.SetLineWidth(2)
        line7.SetLineColor(17)
        line7.Draw()
        line17 = TLine(43.3, ymin, 43.3, ymax)
        line17.SetLineStyle(2)
        line17.SetLineWidth(2)
        line17.SetLineColor(17)
        line17.Draw()
        line8 = TLine(45.5, ymin, 45.5, ymax)
        line8.SetLineStyle(2)
        line8.SetLineWidth(2)
        line8.SetLineColor(17)
        line8.Draw()
        line9 = TLine(50.2, ymin, 50.2, ymax)
        line9.SetLineStyle(2)
        line9.SetLineWidth(2)
        line9.SetLineColor(17)
        line9.Draw()
        line10 = TLine(54.5, ymin, 54.5, ymax)
        line10.SetLineStyle(2)
        line10.SetLineWidth(2)
        line10.SetLineColor(17)
        line10.Draw()
        line11 = TLine(60.65, ymin, 60.65,ymax) #???
        line11.SetLineStyle(2)
        line11.SetLineWidth(2)
        line11.SetLineColor(17)
        line11.Draw()
        line12 = TLine(78.8, ymin, 78.8, ymax) #??
        line12.SetLineStyle(2)
        line12.SetLineWidth(2)
        line12.SetLineColor(17)
        line12.Draw()
        # line13 = TLine(0,1, 100,1) #??
        # line13.SetLineStyle(1)
        # line13.SetLineWidth(2)
        # line13.SetLineColor(13)
        # line13.Draw()
        
        leg = TLegend(0.4,0.2,0.7,0.35);
        #leg.AddEntry(rp1, "\\frac{MC reconstructed}{MC generated true}","LP");
        leg.AddEntry(rp2, "\\frac{MC reconstructed matched to true coordinates}{MC generated true}", "LP")
        leg.SetBorderSize(0);
        leg.SetFillColor(kWhite);
        leg.SetFillStyle(0);
        leg.SetTextSize(0.03);
        leg.Draw("");
        ROOT.SetOwnership(leg,False);

        txt = TPaveText(0.942,0.85,0.942,0.95,"NDC");
        txt.SetFillColor(kWhite);
        txt.SetFillStyle(0);
        txt.SetBorderSize(0);
        txt.SetTextAlign(33);#middle,left
        txt.SetTextFont(42);#helvetica
        txt.SetTextSize(0.03);
        txt.AddText("Simulation this thesis");
        txt.Draw();
        ROOT.SetOwnership(txt,False);

        # txt2 = TPaveText(0.89,0.8,0.89, 0.925,"NDC");
        # txt2.SetFillColor(kWhite);
        # txt2.SetFillStyle(0);
        # txt2.SetBorderSize(0);
        # txt2.SetTextAlign(33);#middle,left
        # txt2.SetTextFont(42);#helvetica
        # txt2.SetTextSize(0.03);
        # txt2.AddText("this thesis");
        # txt2.Draw();
        # ROOT.SetOwnership(txt2,False);

        txt3 = TPaveText(0.93,0.8,0.93,0.925,"NDC");#txt3 = TPaveText(0.92,0.75,0.92,0.90,"NDC");
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


        
            
            
        


if __name__ == "__main__":
    filename = "AnalysisResults_chicall_20240224.root"
    ratiophoton_plot(filename, 1, "20240305/plot_rxy.pdf")
    ratiophoton_plot(filename, 2, "20240305/plot_ratiophoton.pdf")