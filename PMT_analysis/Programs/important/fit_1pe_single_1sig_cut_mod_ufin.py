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
parser.add_option( "-q", "--charge", action="store", type="float", dest="Q_1"  )
parser.add_option( "-s", "--sigma", action="store", type="float", dest="sigma_1"  )
parser.add_option( "-c", "--sigma_cut", action="store", type="string", dest="sigma_number_cuts"  )
(options, args) = parser.parse_args()

path, file = os.path.split(options.filename)
Q_1 = options.Q_1
sigma_1 = options.sigma_1
sigma_n_cuts = options.sigma_number_cuts
basename = os.path.splitext(file)[0]

print "file name: ", basename 
print "path: ", path 
print file 


print "...This file will be used:	", ("%s.root" % basename)
#Get histogram from ROOT file
c1 = ROOT.TCanvas("c1","",800,800)
pad1 = ROOT.TPad("pad1", "pad1", 0, 0.24, 1, 1)
pad2 = ROOT.TPad("pad1", "pad1", 0, 0, 1, 0.22)
pad1.SetBottomMargin(0.00001);
pad1.SetBorderMode(0);
pad1.SetLogy(1);
pad2.SetTopMargin(0.00001);
pad2.SetBottomMargin(0.1);
pad2.SetBorderMode(0);
pad1.Draw("");
pad2.Draw();
pad1.cd();file = ROOT.TFile("%s/%s.root" % (path, basename), "read")
#print "test"
Tree = file.Get("treeData")

print Tree.GetEntries()


hist = ROOT.TH1F("Charge spectrum", "Charge spectrum", 1000,-40,100)
#hist.GetYaxis().SetLabelFont(63);
#hist.GetYaxis().SetLabelSize(16);
#hist.GetYaxis().SetTitle("#");
#hist.Sumw2(1)

for entry in Tree:
	y = entry.peak_integral
	#print y 
	hist.Fill(y)
#print "test"#, ( "%s.root" % file )
hist.Rebin(4)
hist.SetMinimum(1)
#hist.Draw();
#hist.Draw("hist e")

#Set fitting range from maximum of pedestal to last bin with more than 10 counts
left_end = hist.GetXaxis().GetBinCenter(hist.FindFirstBinAbove(1,1))
right_end = hist.GetXaxis().GetBinCenter(hist.FindLastBinAbove(1,1))
#right_end = hist.GetXaxis().GetBinCenter(hist.FindLastBinAbove(500,1))


#Define background function consisting of a gaussian for the pedestal and an exponential for background processes
background = "(((1-[6])/([4]*sqrt(2*3.14159)))*exp(-((x-[3])/(sqrt(2)*[4]))^2)+[6]*(x>[3])*[7]*exp(-[7]*(x-[3])))*exp(-[0])"

#Define background function consisting of a gaussian for the pedestal and an exponential for background processes
background = "(((1-[6])/([4]*sqrt(2*3.14159)))*exp(-((x-[3])/(sqrt(2)*[4]))^2)+[6]*(x>[3])*[7]*exp(-[7]*(x-[3])))*exp(-[0])"

#Define #nr_pe Poisson convoluted Gaussians
gaussian = ""
nr_pe = 15
for ii in range(1,nr_pe):
	gaussian += ("(([0]^%s*exp(-[0]))/%s)*(1/((sqrt(%s)*[2]+[4])*sqrt(2*3.14159)))*exp(-0.5*((x-[1]*%s-[3]-[6]/[7])/(sqrt(%s)*[2]+[4]))^2)+" % (str(ii),str(fac(ii)),str(ii),str(ii),str(ii)))

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
f.SetParameter(1,8.1); # average charge of single pe
f.SetParLimits(1,Q_1,15);
#f.FixParameter(1, 8.9);
f.SetParName(1, "Q_1");
f.SetParameter(2,3.386); # sigma of 1 pe (?)
#~ f.SetParLimits(2,0.5,6.5);
f.SetParLimits(2,2.5,sigma_1);
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
#f.FixParameter(6, 0.05);
f.SetParName(6, "omega")
f.SetParameter(7,0.312); # Coefficient of exponential decrease of thermoelectric proc
f.SetParLimits(7,.001,100000);
#f.SetParLimits(7,2,10);
#f.FixParameter(7, 0.15);
f.SetParName(7, "alpha")

#Fit the response function #nr_fit times to the histogram
nr_fit = 4
for ii in range(nr_fit):
	if(ii==nr_fit-1):
		hist.Fit("f","MERL","");
	else:
		hist.Fit("f","MERL","");
		print "Fit Nr. ", ii, "\r"
		
# get Chisquare and store it in textfile
with open( "%s/chisquare.txt" %path, "w" ) as chifile:
	chi2 = f.GetChisquare() / f.GetNDF()
	print chi2
	chifile.write( "%f" %chi2)

#Extract the fitted paratemeter values to plot the individual components of the PMT response function
f1 = []
for ii in range(1,nr_pe):
	f1.append(ROOT.TF1("f1",("(([0]^%s*exp(-[0]))/%s)*([3]/((sqrt(%s)*[2]+[7])*sqrt(2*3.14159)))*exp(-0.5*((x-[1]*%s-[4]-[5]/[6])/(sqrt(%s)*[2]+[7]))^2)" % (str(ii),str(fac(ii)),str(ii),str(ii),str(ii))),-10, 10000))
	f1[ii-1].SetNpx(10000)
	f1[ii-1].SetLineStyle(3)
	f1[ii-1].SetParameters(f.GetParameter(0),f.GetParameter(1),f.GetParameter(2),f.GetParameter(5),f.GetParameter(3), f.GetParameter(6), f.GetParameter(7), f.GetParameter(4))
	f1[ii-1].SetLineColor(ROOT.kBlue) # Tcolor
	f1[ii-1].Draw("same")

#Plot background contribution
f2 = ROOT.TF1("f2","[3]*(((1-[4])/([2]*sqrt(2*3.14159)))*exp(-((x-[1])/(sqrt(2)*[2]))^2)+[4]*(x>[1])*[5]*exp(-[5]*(x-[1])))*exp(-[0])", -10, 10000);
f2.SetNpx(10000);
f2.SetParameters(f.GetParameter(0),f.GetParameter(3),f.GetParameter(4),f.GetParameter(5),f.GetParameter(6),f.GetParameter(7));
f2.SetLineColor(ROOT.kMagenta);
f2.SetLineStyle(3);
f2.Draw("SAME");

#pad1.SaveAs("%s/pmt_1pe_w_cuts.pdf" % path)
# save Histogram:
with open("%s/histogram_data.txt" % path, "w") as histfile:
	histfile.write("# integral\tcounts")
	for i in range(1, hist.GetXaxis().GetNbins() + 1):
		histfile.write( "%.5f\t%.5f\n" % ( hist.GetBinCenter(i), hist.GetBinContent(i) ) )
		
# save fit function and convolution:
with open("%s/fit_data.txt" % path, "w") as fitfile:
	fitfile.write("# integral\tfit_func\tbackground\tconvolution_1\tconvolution_2\tconvolution_3\tconvolution_4\tconvolution_5\tconvolution_6\tconvolution_7")
	step = 0
	n_steps = 100000
	for i in range(n_steps):
		value = float( -40 + step )
		
		fit_function = f.Eval(value)
		background = f2.Eval(value)
		convolution_1 = f1[0].Eval(value)
		convolution_2 = f1[1].Eval(value)
		convolution_3 = f1[2].Eval(value)
		convolution_4 = f1[3].Eval(value)
		convolution_5 = f1[4].Eval(value)
		convolution_6 = f1[5].Eval(value)
		convolution_7 = f1[6].Eval(value)
		
		fitfile.write( "%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n" % ( value, fit_function, background, convolution_1, convolution_2, convolution_3, convolution_4, convolution_5, convolution_6, convolution_7 ) )
		step += float((100+40))/float(n_steps)


# Calculate residual
pad2.cd();
#x = [ 0 ] * hist.GetXaxis().GetNbins()
#y = [ 0 ] * hist.GetXaxis().GetNbins()
with open( "%s/residual.txt" % path, "w") as g:
	g.write("# integral\tresidual\terror\n")
	h2 = ROOT.TH1F("", "", 250,-40,100)
	h2.GetXaxis().SetLabelFont(63);
	h2.GetXaxis().SetLabelSize(16);
	#h2.GetXaxis().SetTitle("integral [mVns]");
	#h2.GetYaxis().SetTitle("residual [%]");
	#h2.GetYaxis().SetTitleSize(16);
	#h2.GetXaxis().SetTitleSize(16);
	#h2.GetXaxis().SetTitleOffset(-0.9)
	h2.GetYaxis().SetLabelFont(63);
	h2.GetYaxis().SetLabelSize(16);
	h2.GetYaxis().SetRangeUser(-3,3)
	h2.SetStats(0)

	for i in range(hist.GetMaximumBin() - 6, hist.GetXaxis().GetNbins() + 1):
		if ( hist.GetBinContent(i) == 0 ):
			res = 0
			error = 0
			h2.SetBinError( i, 0)
		else:
			#~ res = ( float( hist.GetBinContent(i)) - f.Eval(hist.GetBinCenter(i))) / ( np.sqrt(float( hist.GetBinContent(i))) );     # res of zeroth order poisson according to caio
			res = ( float( hist.GetBinContent(i)) - f.Eval(hist.GetBinCenter(i))) / ( f.Eval(hist.GetBinCenter(i)) );
		#h2.SetBinError( i, np.sqrt( ( np.sqrt(hist.GetBinContent(i))/f.Eval(hist.GetBinCenter(i)) )**2 + ( ( (hist.GetBinContent(i)) * np.sqrt( f.Eval(hist.GetBinCenter(i)) ) )/( (f.Eval(hist.GetBinCenter(i)))**2 ) )**2 ));
			#~ h2.SetBinError( i, np.sqrt(  ( ( np.sqrt( hist.GetBinContent(i) ) - ( hist.GetBinContent(i) - f.Eval(hist.GetBinCenter(i)) ) * 0.5 * 1/np.sqrt(hist.GetBinContent(i)) ) / (np.sqrt(hist.GetBinContent(i))) )**2 ))		# error of zeroth order poisson
			h2.SetBinError( i, np.sqrt(hist.GetBinContent(i)) / f.Eval(hist.GetBinCenter(i)) )
		#error  = np.sqrt( ( np.sqrt(hist.GetBinContent(i))/f.Eval(hist.GetBinCenter(i)) )**2 + ( ( (hist.GetBinContent(i)) * np.sqrt( f.Eval(hist.GetBinCenter(i)) ) )/( (f.Eval(hist.GetBinCenter(i)))**2 ) )**2 )
			error = np.sqrt(hist.GetBinContent(i)) / f.Eval(hist.GetBinCenter(i))
			#~ error = np.sqrt(  ( ( np.sqrt( hist.GetBinContent(i) ) - ( hist.GetBinContent(i) - f.Eval(hist.GetBinCenter(i)) ) * 0.5 * 1/np.sqrt(hist.GetBinContent(i)) ) / (np.sqrt(hist.GetBinContent(i))) )**2 )		# error of zeroth order poisson
		#print hist.GetBinContent(i), "\t", f.Eval(hist.GetBinCenter(i))
			#print ( float( hist.GetBinContent(i)) - f.Eval(hist.GetBinCenter(i))) / float(hist.GetBinContent(i));
		#x[i-1] = float(hist.GetBinCenter(i))
		#y[i-1] = res
			g.write( "%.5f\t%.5f\t%.5f\n" % ( hist.GetBinCenter(i), res, error ) )
		h2.SetBinContent(i, res);
	g.write("\n")
	h2.Draw();
	c1.cd();
#x = np.array(x)
#y = np.array(y)
#graph = ROOT.TGraph(len(x))
#c2 = ROOT.TCanvas("c2","",2500,2200)
#graph = ROOT.TGraph(len(x), x, y)
#graph.Draw();
#c2.SetLogy(0)
#c2.SaveAs("%s/residual.pdf" % path)

c1.SaveAs("%s/pmt_1pe_res_1sig_cut_unmod_final.pdf" % path)



# write fit parameters:
with open( "%s/parameters.txt" % (path), "w") as g:
	#~ g.write("#mu\tmu_error\tQ1\tQ1_error\tsigma_1\tsigma_1_error\tQ0\tQ0_error\tsigma_0\tsigma_0_error\tscaling\tscaling_error\tw\tw_error\talpha\talpha_error\n")
	g.write( str( hist.GetEntries() ) )
	g.write( "\t" )
	g.write( str( chi2 ) )
	g.write( "\t" )
	for i in range(8):
		g.write(str(f.GetParameter(i)))
		g.write("\t")
		g.write(str(f.GetParError(i)))
		g.write("\t")

