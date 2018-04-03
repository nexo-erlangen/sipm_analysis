import sys
import os
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
from difflib import Differ
#~ ROOT.gROOT.SetBatch(True)
#~ ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")
#~ ROOT.gStyle.SetOptFit(111);
#~ fac = ROOT.TMath.Factorial
parser = OptionParser()

path = "/home/vault/capm/sn0527/PMT_LED_Temp/Auswertung_diff_cuts/"
temp_list = [ "+23", "-9", "-43", "-73" ]
file_list = [ "0_sigma", "1_sigma", "2_sigma", "3_sigma", "no_cut" ]
Temp = [23.0, -9.0, -43.0, -73.0]
w_mean_list_all = []

# reading textfile with parameters and calculating weighted mean with error
# loop over all sigma
for sig in range(len(file_list)):
	#~ for s in range(2):
	w_mean_list = []
	# loop over different temperatures
	for i in range(len(temp_list)):
		parameters = np.array( np.genfromtxt( "%sparameter_T=%sC_%s.txt" % (path, temp_list[i], file_list[sig] ), delimiter = "\t" ) )
		#~ parameters = np.loadtxt( "%sparameter_T=%sC_%s.txt" % (path, temp_list[i], file_list[int(s)] ), skiprows=1 )
		mu = parameters[:, 2]
		mu_error = parameters[:, 3]
		Q1 = parameters[:, 4]
		Q1_error = parameters[:, 5]

		Q1_std_weighted = np.sqrt(1./np.nansum([1./s**2 for s in Q1_error]))
		Q1_weighted = ( np.nansum([Q1[n]/(Q1_error[n]**2) for n in range(len(Q1))]) )/( np.nansum([1./(Q1_error[n]**2) for n in range(len(Q1))]) )
		line = []
		line.append( "%.5f" % Temp[i] )
		line.append( "%.5f" % Q1_weighted )
		line.append( "%.5f" % Q1_std_weighted )
		line.append( "%.5f" % ( Q1_std_weighted/Q1_weighted ) )
		#~ print line
		w_mean_list.append( line )
	w_mean_list_all.append(w_mean_list)
#~ print w_mean_list_all
w_mean_list_all = np.array(w_mean_list_all, dtype=float)
#~ print w_mean_list_all

'''
# plot weighted mean with error in dependence of temperature for different sigma
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("%s" % "comparison different cuts" )    
ax1.set_xlabel('T [C]')
ax1.set_ylabel('%s' % "weighted mean Q1" )
for s in range(len(file_list)):
	ax1.errorbar( w_mean_list_all[s][:,0],w_mean_list_all[s][:,1], yerr=w_mean_list_all[s][:,2], ls="none", marker="o", label='%s' % file_list[s] )
leg = ax1.legend()
fig.show()
raw_input("")
'''


# fitting curves with ROOT
# set no display of canvases 
ROOT.gROOT.SetBatch(True)

chi2 = []
NDF = []
p0 = []
p0_err = []
p1 = []
p1_err = []

# loop over different sigma
for s in range(len(file_list)):
#~ for s in range(2):
	x = w_mean_list_all[s][:,0].tolist()
	y = w_mean_list_all[s][:,1].tolist()
	xerr = [0] * 15
	yerr = w_mean_list_all[s][:,3].tolist()

	c1 = ROOT.TCanvas("c1","",800,800)
	graph1 = ROOT.TGraphErrors( len(x), np.array(x), np.array(y), np.array(xerr), np.array(yerr) )

	graph1.Draw()
	#~ raw_input("")

	# f(x) = p0 + p1*x
	function = "pol1"
	f = ROOT.TF1("f",function)
	graph1.Fit(f)
	f.Draw("same")

	# get parameters
	chisquare = f.GetChisquare()
	ndof = f.GetNDF()
	param0 = f.GetParameter(0)
	param0_err = f.GetParError(0)
	param1 = f.GetParameter(1)
	param01_err = f.GetParError(1)
	
	chi2.append( f.GetChisquare() )
	NDF.append( f.GetNDF() )
	p0.append( f.GetParameter(0) )
	p0_err.append( f.GetParError(0) )
	p1.append( f.GetParameter(1) )
	p1_err.append( f.GetParError(1) )
	
print chi2


# plot weighted mean with error in dependence of temperature for different sigma
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("%s" % "comparison different cuts" )    
ax1.set_xlabel('T [C]')
ax1.set_ylabel('%s' % "weighted mean Q1" )
xvalue = np.arange( -120, 40, 0.01 )
Q1_m100 = []
rel_var_from = []
for s in range(len(file_list)):
	ax1.errorbar( w_mean_list_all[s][:,0],w_mean_list_all[s][:,1], yerr=w_mean_list_all[s][:,2], ls="none", marker="o", label='%s' % file_list[s] )
	yvalue = p0[s] + p1[s] * xvalue
	plt.plot( xvalue, yvalue )
	# extrapolate each curve to -100C
	Q1_m100.append(p0[s] + p1[s] * (-100.0))
	#~ print Q1_m100[s]
	# relative variation from 0sigma
	rel_var_from.append( ( ( p0[0] + p1[0] * (-100.0) ) - ( p0[s] + p1[s] * (-100.0) ) ) / ( p0[0] + p1[0] * (-100.0) ) )
	#~ print rel_var_from[s]
leg = ax1.legend()
#~ fig.show()
print Q1_m100

# plot variation of Q1 depending of cut
fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.set_title("%s" % "comparison different cuts" )    
ax1.set_xlabel('cut')
ax1.set_ylabel('%s' % "weighted mean Q1(-100C)" )
for s in range(len(file_list)):
	ax1.plot( s, Q1_m100[s], ls="none", marker="o", label="none" )
#~ leg = ax1.legend()
ax2 = fig.add_subplot(212)
#~ ax2.set_title("%s" % "deviation from 0sigma cut" )    
ax2.set_xlabel('cut')
ax2.set_ylabel('%s' % "rel deviation" )
for s in range(len(file_list)):
	ax2.plot( s, rel_var_from[s], ls="none", marker="o", label='%s_sigma' % s )
	
#~ fig.show()

# calculate PDE with different cut results for Q1
A_pmt = 420.25 #mm^2
A_sipm = 36.0 # mm^2
#DE = 0.29 # sagt Patrick # 0.7*0.311=0.22
DE = 0.2247
R = 0.8592 # fuer Al 60-90, Cu 20-30 reflectivity


delta_DE = 0.03
delta_A_pmt = 0.4
delta_A_sipm = 0.2
delta_R = 0.00624
delta_1pe_PMT = 0.106602126675181 # hoechstwahrscheinlich

path_2 = "/home/vault/capm/sn0527/FBK-LF-STD_Pos1/Auswertung_PDE/"

integrals = np.array( np.genfromtxt( "%sintegrals_for_pde.txt" % (path_2 ), delimiter = "\t" ) )
print integrals

fig_pde = plt.figure()
ax_pde = fig_pde.add_subplot(111)
#~ ax_pde.set_ylim( [0,0.12] )
for s in range(len(file_list)):
	ax_pde.plot( ( integrals[:,0] + 29.558)*(-1), ( integrals[:,1] / integrals[:,3] ) / ( integrals[:,5] / Q1_m100[s] ) *A_pmt / A_sipm * DE / R, label='%s_sigma' % s )
leg = ax_pde.legend()
fig_pde.show()

raw_input("")

