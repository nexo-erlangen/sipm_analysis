import sys
import os
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser
from difflib import Differ
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")
ROOT.gStyle.SetOptFit(111);
fac = ROOT.TMath.Factorial
parser = OptionParser()

param_list = [ "Q1", "sigma1", "Q0", "sigma0", "scaling", "omega", "alpha" ]
temp_list = [ "+23", "-9", "-43", "-73" ]
file_list = [ "0_sigma", "1_sigma", "2_sigma", "3_sigma", "no_cut" ]
c_list = [ "# entries", "reduced chi2" ]
sigma_list = []


# reading
# loop over sigmas
for j in range(5):
	#~ print j
	data_temperatures = []
	# loop over temperatures
	for i in range(4):
		data = np.loadtxt( "/home/vault/capm/sn0527/PMT_LED_Temp/Auswertung_diff_cuts/parameter_T=%sC_%s.txt" % (temp_list[i], file_list[j] ), skiprows=1 )
		data_temperatures.append(data)
	sigma_list.append(data_temperatures)

'''
# test for reading data
# loop over temperatures
for i in range(4):
	data = np.loadtxt( "/home/vault/capm/sn0527/PMT_LED_Temp/Auswertung_diff_cuts/parameter_T=%sC_3_sigma.txt" % temp_list[i], skiprows=1 )
	data_temperatures.append(data)
	
data_temperatures = np.array(data_temperatures)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_title("%s sigma" % str(1) )    
ax1.set_xlabel('mu (average number of detected photons)')
ax1.set_ylabel('%s' % "test")
ax1.plot( data_temperatures[0][:,2],data_temperatures[0][:,1], ls="none", marker="o", label='%s sigma' % str(1) )
leg = ax1.legend()
fig.show()
raw_input("")
'''

# plot all parameter in dependence of mu for different temperatures in one plot, for every sigma one plot
data_temperatures = np.array(data_temperatures)

# loop over different sigmas
for s in range(5):
	# loop over different parameters
	n = 0
	for p in range(7):
		fig = plt.figure()
		ax1 = fig.add_subplot(111)
		ax1.set_title("%s" % file_list[s] )    
		ax1.set_xlabel('mu (average number of detected photons)')
		ax1.set_ylabel('%s' % param_list[p] )
		# loop over different temperatures
		for i in range(4):
			ax1.errorbar( sigma_list[s][i][:,2],sigma_list[s][i][:,4 + n ], xerr=sigma_list[s][i][:,3], yerr=sigma_list[s][i][:,5 + n ], ls="none", marker="o", label='T = %sC' % temp_list[i] )
		n += 2
		leg = ax1.legend()
		fig.show()
raw_input("")


'''
# plot chi2 and entries in dependence of mu for different sigma in one plot, for every temperature one plot
# loop over parameters
for c in range(2):
	# loop over all temperatures
	for t in range(4):
		fig = plt.figure()
		ax1 = fig.add_subplot(111)
		ax1.set_title("%s" % temp_list[t] )    
		ax1.set_xlabel('mu (average number of detected photons)')
		ax1.set_ylabel('%s' % c_list[c] )
		# loop over different sigmas
		for s in range(5):
			ax1.plot( sigma_list[s][t][:,2],sigma_list[s][t][:,c], ls="none", marker="o", label='%s' % file_list[s] )
			leg = ax1.legend()
			fig.show()
raw_input("")
'''
'''
# plot all parameters in dependence of mu for different sigma in one plot, for every temperature one plot
# loop over all temperatures
for t in range(4):
	# loop over parameters
	n = 0
	for p in range(7):
		fig = plt.figure()
		ax1 = fig.add_subplot(111)
		ax1.set_title("%s" % temp_list[t] )    
		ax1.set_xlabel('mu (average number of detected photons)')
		ax1.set_ylabel('%s' % param_list[p] )
		#~ ax1.errorbar( sigma_list[s][t][:,2],sigma_list[s][t][:,4 + n ], xerr=sigma_list[s][t][:,3], yerr=sigma_list[s][t][:,5 + n ], ls="none", marker="o", label='%s sigma' % str(s) )
		# loop over different sigmas
		for s in range(5):
			ax1.errorbar( sigma_list[s][t][:,2],sigma_list[s][t][:,4 + n ], xerr=sigma_list[s][t][:,3], yerr=sigma_list[s][t][:,5 + n ], ls="none", marker="o", label='%s' % file_list[s] )
		n += 2
		leg = ax1.legend()
		fig.show()
raw_input("")
'''

exit()