import sys
import os
import ROOT
import numpy as np
import matplotlib as plot
from optparse import OptionParser
from difflib import Differ
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")
ROOT.gStyle.SetOptFit(111);
fac = ROOT.TMath.Factorial
parser = OptionParser()

#When calling this script the (path/)filename should be passed as option
parser.add_option("-f", "--file",action="store", type="string", dest="filename")
(options, args) = parser.parse_args()

path, file = os.path.split(options.filename)
basename = os.path.splitext(file)[0]

print "file name: ", basename 
print "path: ", path 
print file 


print "...This file will be used:	", ("%s.root" % basename)
#Get histogram from ROOT file
c1 = ROOT.TCanvas("c1","",2500,2200)
file = ROOT.TFile("%s/%s.root" % (path, basename), "read")
#print "test"
Tree = file.Get("treeData")

print Tree.GetEntries()


hist = ROOT.TH1F("hist", "hist", 1000,-1,1)
#hist.Sumw2(1)

for entry in Tree:
	y = entry.baseline_rms
	#print y 
	hist.Fill(y)
#print "test"#, ( "%s.root" % file )
hist.Rebin(4)
hist.SetMinimum(1)
hist.Draw();
#hist.Draw("hist e")

#Save plot as pdf
c1.SetLogy(1)
c1.SaveAs("%s/pmt_1pe_datatree_histogram_baseline_rms.pdf" % path)