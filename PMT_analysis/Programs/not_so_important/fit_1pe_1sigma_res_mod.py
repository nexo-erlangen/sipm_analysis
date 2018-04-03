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

def main():
	# zu verwendete Runs:
	#run_list=["00703", "00704", "00705"]
	run_list = [ [ "00695", "00696", "00697", "00698"], ["00700", "00701", "00702"], ["00703", "00704", "00705"], ["00706", "00707", "00708", "00709"] ]
	Temp_list = [ "+23", "-73", "-43", "-9"]
	#~ Q1_list = [ [ 0.5, 0.5, 7 ], [ 7, 0.5, 0.5 ], [ 0.5, 0.5, 0.5 ], [ 0.5, 0.5, 0.5, 0.5 ] ]
	u_sigma_list = [ [ 2.9, 6.5, 6.5, 3.5 ], [ 3.4, 3.7, 6.5 ], [ 6.5, 3.7, 6.5 ], [ 3.4, 3.5, 6.5, 6.5 ] ]

		
	# loop over temperatures
	#~ for j in range(4):
	for j in range(4):
		print "\n"
		print "--------------------------- new temperature ---------------------------"
		with open( "/home/vault/capm/sn0527/PMT_LED_Temp/Auswertung_mod/limited/parameter_1sig_lim_T=%sC.txt" % (Temp_list[j]  ), "w" ) as plotfile:
			plotfile.write( "#entries\tchi2\tmu\tmu_error\tQ1\tQ1_error\tsigma_1\tsigma_1_error\tQ0\tQ0_error\tsigma_0\tsigma_0_error\tscaling\tscaling_error\tw\tw_error\talpha\talpha_error\n" )
			# loop over intensities
			for i in range(len(run_list[j])):
				print "\n"
				print "--------------------------- new Run ---------------------------"
				print "\n"
				parameters = ""
				# fit with limited parameters:
				os.system( "python /home/vault/capm/sn0527/PMT_LED_Temp/fit_1pe_single_res_mod.py --file /home/vault/capm/sn0527/PMT_LED_Temp/Run%s/fit_tree_1_sigma.root -s %f" %(  run_list[j][i], u_sigma_list[j][i] )  )				
				os.system( "cp /home/vault/capm/sn0527/PMT_LED_Temp/Run%s/pmt_1pe_w_cuts_res_mod_temp.pdf /home/vault/capm/sn0527/PMT_LED_Temp/Auswertung_mod/limited/pmt_1pe_1sig_lim_Run%s.pdf" %( run_list[j][i], run_list[j][i] ) )
				print "--------------------------- copy pdf ---------------------------"
				with open( "/home/vault/capm/sn0527/PMT_LED_Temp/Run%s/parameters.txt" %run_list[j][i], "r") as paramfile:
					parameters = paramfile.readline()
				plotfile.write( "%s" % parameters )
				print "--------------------------- writing parameters ---------------------------"
				plotfile.write( "\n" )
				
				# try to fit with free parameters:
				#~ os.system( "python /home/vault/capm/sn0527/PMT_LED_Temp/fit_1pe_single_res_mod.py --file /home/vault/capm/sn0527/PMT_LED_Temp/Run%s/fit_tree_1_sigma.root" %(  run_list[j][i] )  )				
				#~ os.system( "cp /home/vault/capm/sn0527/PMT_LED_Temp/Run%s/pmt_1pe_w_cuts_res_mod_temp.pdf /home/vault/capm/sn0527/PMT_LED_Temp/Auswertung_mod/free_fit/pmt_1pe_1sig_free_Run%s.pdf" %( run_list[j][i], run_list[j][i] ) )
				#~ print "--------------------------- copy pdf ---------------------------"
				#~ with open( "/home/vault/capm/sn0527/PMT_LED_Temp/Run%s/parameters.txt" %run_list[j][i], "r") as paramfile:
					#~ parameters = paramfile.readline()
				#~ plotfile.write( "%s" % parameters )
				#~ print "--------------------------- writing parameters ---------------------------"
				#~ plotfile.write( "\n" )
				
							
	
if __name__== "__main__":
	main()