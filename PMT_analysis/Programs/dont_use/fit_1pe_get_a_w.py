import sys
import os
import ROOT
import numpy as np
import matplotlib as plot
import time
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
	start_time = time.clock()
	# zu verwendete Runs:
	result_list=[]
	list_all_error=[]
	run_list=["00695", "00696", "00697", "00698"]
	omega_left = 0.025
	omega_right = 0.07
	steps = 30
	#omega_left = 0.01
	#omega_right = 0.05
	#steps = 8
	
	alpha_left = 0.1
	alpha_right = 0.55
	alpha_steps = 30
	#alpha_left = 0.1
	#alpha_right = 0.4
	#alpha_steps = 8
	#counter = 0
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
		# Variation von alpha:
		param_list=[]
		for alpha in np.linspace(alpha_left, alpha_right, alpha_steps):
			# Schleife zur Variation von omega:
			#omega_list=[]
			#print "\n alpha = ", alpha, "\n"
			for omega in np.linspace(omega_left, omega_right, steps):
				#counter +=  1
				#print "\nStep:\t", counter, "from\t", steps*alpha_steps*len(run_list)
				#print "<<<<< ", float(counter)/(steps*alpha_steps*len(run_list)), " % <<<<<"
				#print "\n"
				fit_list = fit_histogram(hist, omega, alpha, path)
				if fit_list[-1] > 20.0:
					param_list.append([alpha, omega, np.nan, np.nan, np.nan]) 
					#continue # 0er reinschreiben verfaelscht/ vergroessert Mittelwert
				else:
					param_list.append([alpha, omega, fit_list[2], fit_list[3], fit_list[-1] ]) 
			#Q1_list = np.array(omega_list)[:,2]
			#Q1err_list = np.array(omega_list)[:,3]
			#alpha_list.append(Q1_list)
		result_list.append(param_list)	
		#list_all_error.append(alpha_list)
		file.Close()
			
	with open("a_o_Q1_T=+23C_013.dat", "w") as f:
		f.write("# alpha\tomega\tmean(Q1)\tmeanw(Q1)\tstd(Q1)\tstdw(Q1)\n")
		for j in range(alpha_steps):
			for i in range(steps):
				alpha, omega = result_list[0][j*steps + i][0:2]
				Q1=np.array(result_list)[:,j*steps + i,2]
				Q1_mean = np.nanmean( Q1 )
				Q1_std = np.nanstd( Q1 )
				Q1_sigma = np.array(result_list)[:,j*steps + i,3]
				Q1_std_weighted = np.sqrt(1./np.nansum([1./s**2 for s in Q1_sigma]))
				Q1_weighted = ( np.nansum([Q1[n]/(Q1_sigma[n]**2) for n in range(len(Q1))]) )/( np.nansum([1./(Q1_sigma[n]**2) for n in range(len(Q1))]) )
				#sum = 0
				#for m in range(len(run_list)):
					#sum += 1.0/(Q1_simga[m])**2
				#Q1_std_weighted = np.sqrt(1.0/sum)
				# chi_red = 
				# if chi_red < 55.0:
				if float(Q1_std)>0.00001:
					f.write("%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n" %(alpha, omega,Q1_mean, Q1_weighted , Q1_std, Q1_std_weighted))
			f.write("\n")
		f.write("\n")
	stop_time = time.clock()
	print "\n\n elapsed time:\t", stop_time - start_time, " s"
	print "\n time per fit:\t", (stop_time - start_time)/(alpha_steps*steps), " s"
	

def fit_histogram(hist, val_omega, val_alpha, path):
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
	f.SetParLimits(1,4,15);
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
	#f.SetParLimits(7,.001,100000);
	f.FixParameter(7, val_alpha);
	f.SetParName(7, "alpha")

	#Fit the response function #nr_fit times to the histogram
	nr_fit = 10
	for ii in range(nr_fit):
		if(ii==nr_fit-1):
			hist.Fit("f","MERL","");
		else:
			hist.Fit("f","MERL","");
			print "Fit Nr. ", ii, "\r"
			print "Run: ", path

	#Extract the fitted paratemeter values to plot the individual components of the PMT response function
	f1 = []
	for ii in range(1,nr_pe):
		f1.append(ROOT.TF1("f1",("(([0]^%s*exp(-[0]))/%s)*([3]/([2]*sqrt(2*%s*3.14159)))*exp(-((x-[1]*%s-[4]-[5]/[6])/(sqrt(2*%s)*[2]))^2)" % (str(ii),str(fac(ii)),str(ii),str(ii),str(ii))),-10, 10000))
		f1[ii-1].SetNpx(10000)
		f1[ii-1].SetLineStyle(3)
		f1[ii-1].SetParameters(f.GetParameter(0),f.GetParameter(1),f.GetParameter(2),f.GetParameter(5),f.GetParameter(3), f.GetParameter(6), f.GetParameter(7))
		f1[ii-1].SetLineColor(ROOT.kBlue-4-ii) # Tcolor
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
	
	fitparameters.append(f.GetChisquare()/f.GetNDF()) # soll kleiner sein als 10 zB
	return fitparameters
	#liste= np.array(f.GetParameters())
	#print liste
	
if __name__ == "__main__":
	main()
	