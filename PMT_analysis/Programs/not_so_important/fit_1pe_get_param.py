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
'''parser.add_option("-f", "--file",action="store", type="string", dest="filename")
(options, args) = parser.parse_args()

path, file = os.path.split(options.filename)
basename = os.path.splitext(file)[0]
'''

def main():
	# zu verwendete Runs:
	run_list=["00703", "00704", "00705"]
	result_list=[]
	for run in run_list:
		result_list.append(fit_histogram(run))
	
	with open("fit_parameters_T=-43C.dat", "w") as f:
		f.write("#mu\tmu_error\tQ1\tQ1_error\tsigma_1\tsigma_1_error\tQ0\tQ0_error\tsigma_0\tsigma_0_error\tscaling\tscaling_error\tw\tw_error\talpha\talpha_error\n")
		for result in result_list:
			for entry in result:
				f.write("%.5f\t" % entry)
			f.write("\n")
		f.write("\n")
		

def fit_histogram(run_number):
	path = "/home/vault/capm/sn0527/PMT_LED_Temp/Run" + run_number + "/Auswertung/integral_15-35ns"
	file = "hist.root"

	print "path: ", path 
	print file 


	print "...This file will be used:	", file
	#Get histogram from ROOT file
	c1 = ROOT.TCanvas("c1","",2500,2200)
	file = ROOT.TFile("%s/%s" % (path, file), "read") # hier passiert ansch der Fehler
	#print "test"
	Tree = file.Get("histogramm")

	print Tree.GetEntries()


	hist = ROOT.TH1F("hist", "hist", 1000,-20,70)
	#hist.Sumw2(1)

	for entry in Tree:
		y = entry.integral
		#print y 
		hist.Fill(y)
	#print "test"#, ( "%s.root" % file )
	hist.Rebin(4)
	hist.SetMinimum(1)
	hist.Draw();
	#hist.Draw("hist e")

	#Set fitting range from maximum of pedestal to last bin with more than 10 counts
	left_end = hist.GetXaxis().GetBinCenter(hist.GetMaximumBin() - 9)
	right_end = hist.GetXaxis().GetBinCenter(hist.FindLastBinAbove(15,1))
	#right_end = hist.GetXaxis().GetBinCenter(hist.FindLastBinAbove(500,1))


	#Define background function consisting of a gaussian for the pedestal and an exponential for background processes
	background = "(((1-[6])/([4]*sqrt(2*3.14159)))*exp(-((x-[3])/(sqrt(2)*[4]))^2)+[6]*(x>[3])*[7]*exp(-[7]*(x-[3])))*exp(-[0])"

	#Define #nr_pe Poisson convoluted Gaussians
	gaussian = ""
	nr_pe = 15
	for ii in range(1,nr_pe):
		gaussian += ("(([0]^%s*exp(-[0]))/%s)*(1/([2]*sqrt(2*%s*3.14159)))*exp(-((x-[1]*%s-[3]-[6]/[7])/(sqrt(2*%s)*[2]))^2)+" % (str(ii),str(fac(ii)),str(ii),str(ii),str(ii)))

	gaussian = gaussian[:-1]

	#Combine both components of the PMT response function to one function
	function = "[5]*(" + background + "+" + gaussian +")"


	f = ROOT.TF1("f",function,left_end,right_end)
	f.SetLineColor(ROOT.kRed);
	f.SetNpx(100000);

	#Define parameter ranges and start values
	f.SetParameter(0,0.1129); # number of primary pe
	f.SetParLimits(0,0,10);
	#f.FixParameter(0, 0.03805);
	f.SetParName(0, "mu");
	f.SetParameter(1,8.931); # average charge of single pe
	f.SetParLimits(1,7,15);
	#f.FixParameter(1, 8.9);
	f.SetParName(1, "Q_1");
	f.SetParameter(2,3.386); # sigma of 1 pe (?)
	f.SetParLimits(2,0.5,6.5);
	#f.FixParameter(2,3.386f);
	f.SetParName(2,"sigma_1")
	#~ f.SetParameter(3,hist.GetXaxis().GetBinCenter(hist.GetMaximumBin()));
	#~ f.SetParLimits(3,hist.GetXaxis().GetBinCenter(hist.GetMaximumBin())*.05,hist.GetXaxis().GetBinCenter(hist.GetMaximumBin())*3.5);
	f.SetParameter(3, 0.007333); # Mittelwert von Pedestal
	f.SetParLimits(3, -2, 12);
	#f.FixParameter(3, 0.007279);
	f.SetParName(3, "Q_0")
	f.SetParameter(4,0.9019); # standard deviation of pedestal
	f.SetParLimits(4,0,7.7);
	#f.FixParameter(4,0.9019);
	f.SetParName(4, "sigma_0")
	f.SetParameter(5,2.48E5); # Skalierung von Spektrum
	f.SetParLimits(5, 0, 10000000);
	#f.FixParameter(5, 2.48E5);
	f.SetParName(5, "Scaling")
	f.SetParameter(6,5.316E-8); # Probability of accompanied thermal electron emission process
	f.SetParLimits(6,0,1);
	#f.FixParameter(6, 0.025);
	f.SetParName(6, "omega")
	f.SetParameter(7,0.312); # Coefficient of exponential decrease of thermoelectric proc
	f.SetParLimits(7,.001,100000);
	#f.FixParameter(7, 0.312);
	f.SetParName(7, "alpha")

	#Fit the response function #nr_fit times to the histogram
	nr_fit = 20
	for ii in range(nr_fit):
		if(ii==nr_fit-1):
			hist.Fit("f","MERL","");
		else:
			hist.Fit("f","MERL","");
			print "Fit Nr. ", ii, "\r"

	#Extract the fitted paratemeter values to plot the individual components of the PMT response function
	f1 = []
	for ii in range(1,nr_pe):
		f1.append(ROOT.TF1("f1",("(([0]^%s*exp(-[0]))/%s)*([3]/([2]*sqrt(2*%s*3.14159)))*exp(-((x-[1]*%s-[4]-[5]/[6])/(sqrt(2*%s)*[2]))^2)" % (str(ii),str(fac(ii)),str(ii),str(ii),str(ii))),-10, 10000))
		f1[ii-1].SetNpx(10000)
		f1[ii-1].SetLineStyle(3)
		f1[ii-1].SetParameters(f.GetParameter(0),f.GetParameter(1),f.GetParameter(2),f.GetParameter(5),f.GetParameter(3), f.GetParameter(6), f.GetParameter(7))
		f1[ii-1].SetLineColor(ROOT.kBlue+5-ii) # Tcolor
		f1[ii-1].Draw("same")

	#Plot background contribution
	f2 = ROOT.TF1("f2","[3]*(((1-[4])/([2]*sqrt(2*3.14159)))*exp(-((x-[1])/(sqrt(2)*[2]))^2)+[4]*(x>[1])*[5]*exp(-[5]*(x-[1])))*exp(-[0])", -10, 10000);
	f2.SetNpx(10000);
	f2.SetParameters(f.GetParameter(0),f.GetParameter(3),f.GetParameter(4),f.GetParameter(5),f.GetParameter(6),f.GetParameter(7));
	f2.SetLineColor(ROOT.kMagenta);
	f2.SetLineStyle(3);
	f2.Draw("SAME");

	#Save plot as pdf
	c1.SetLogy(1)
	c1.SaveAs("%s/pmt_1pe.pdf" % path)

	# Get fit parameters
	fitparameters=[]
	for i in range(8):
		fitparameters.append(f.GetParameter(i))
		fitparameters.append(f.GetParError(i))

	return fitparameters
	#liste= np.array(f.GetParameters())
	#print liste
	f.Close()
	
if __name__ == "__main__":
	main()
	