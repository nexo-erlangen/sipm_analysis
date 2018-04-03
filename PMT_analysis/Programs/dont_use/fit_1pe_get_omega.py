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
	list_all_Q1=[]
	list_all_error=[]
	run_list=["00695", "00696", "00697", "00698"]
	omega_left = 0.024
	omega_right = 0.028
	steps = 50
	for run in run_list:
		path = "Run" + run + "/Auswertung/integral_15-35ns"
		file = "hist.root"

		print "path: ", path 
		print file 


		print "...This file will be used:	", file
		#Get histogram from ROOT file
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
		# Schleife zur Variation von omega:
		result_list=[]
		for omega in np.linspace(omega_left, omega_right, steps):
			result_list.append(fit_histogram(hist, omega, path))
		Q1_list = np.array(result_list)[:,2]
		Q1err_list = np.array(result_list)[:,3]
		list_all_Q1.append(Q1_list)	
		list_all_error.append(Q1err_list)
		file.Close()
			
	with open("omega_Q1_T=+23C.dat", "w") as f:
		f.write("# omega\tmean(Q1)\tstd(Q1)\n")	
		for i, omega in enumerate(np.linspace(omega_left, omega_right, steps)):
			Q1 = np.array(list_all_Q1)[:,i]
			f.write("%.5f\t%.5f\t%.5f\n" %(omega, np.mean( Q1 ), np.std( Q1 ) ) )
		f.write("\n")					

	

def fit_histogram(hist, val_omega, path):
	c1 = ROOT.TCanvas("c1","",2500,2200)
	
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
	#f.SetParLimits(6,0,1);
	f.FixParameter(6, val_omega);
	f.SetParName(6, "omega")
	f.SetParameter(7,0.312); # Coefficient of exponential decrease of thermoelectric proc
	f.SetParLimits(7,.001,100000);
	#f.FixParameter(7, 0.312);
	f.SetParName(7, "alpha")

	#Fit the response function #nr_fit times to the histogram
	nr_fit = 10
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
	
if __name__ == "__main__":
	main()
	