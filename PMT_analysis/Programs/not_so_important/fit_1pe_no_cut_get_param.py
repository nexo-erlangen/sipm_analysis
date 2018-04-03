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
	run_list = [ ["00695", "00696", "00697", "00698"], ["00700", "00701", "00702"], ["00703", "00704", "00705"], ["00706", "00707", "00708", "00709"] ]
	#sigma1_limit_list = [5.5, 5.2, 4.8]
	sigma1_limit_list = [ [ 4.2, 4.3, 3.9, 4.6 ], [ 5.0, 5.3, 5.2 ], [5.0, 5.2, 4.8], [ 4.9, 4.8, 4.6, 4.5 ] ]
	Q1_value = [ 7.0, 9.0, 9.0, 7.0 ]
	Temp_list = [ "+23", "-73", "-43", "-9"]

		
	# loop over temperatures
	#~ for j in range(4):
	for j in range(4):
		print "\n"
		print "--------------------------- new temperature ---------------------------"
		with open( "/home/vault/capm/sn0527/PMT_LED_Temp/Auswertung_diff_cuts/parameter_T=%sC_no_cut.txt" % (Temp_list[j]  ), "w" ) as plotfile:
			plotfile.write( "#entries\tchi2\tmu\tmu_error\tQ1\tQ1_error\tsigma_1\tsigma_1_error\tQ0\tQ0_error\tsigma_0\tsigma_0_error\tscaling\tscaling_error\tw\tw_error\talpha\talpha_error\n" )
			# loop over intensities
			for i in range(len(run_list[j])):
				print "\n"
				print "--------------------------- new Run ---------------------------"
				print "\n"
				parameters = ""
				os.system( "python /home/vault/capm/sn0527/PMT_LED_Temp/fit_1pe_single_res_no_cut.py --file /home/vault/capm/sn0527/PMT_LED_Temp/Run%s/datatree.root -q %f -s %f -c %s" %(  run_list[j][i], Q1_value[j], sigma1_limit_list[j][i], "no_cut" )  )
				os.system( "cp /home/vault/capm/sn0527/PMT_LED_Temp/Run%s/pmt_1pe_test_sigma_res_diff_cuts.pdf /home/vault/capm/sn0527/PMT_LED_Temp/Auswertung_diff_cuts/pmt_1pe_Run%s_no_cut.pdf" %( run_list[j][i], run_list[j][i] ) )
				print "--------------------------- copy pdf ---------------------------"
				with open( "/home/vault/capm/sn0527/PMT_LED_Temp/Run%s/parameters.txt" %run_list[j][i], "r") as paramfile:
					parameters = paramfile.readline()
				plotfile.write( "%s" % parameters )
				print "--------------------------- writing parameters ---------------------------"
				plotfile.write( "\n" )
				
							
#	os.system( "python /home/vault/capm/sn0527/PMT_LED_Temp/fit_1pe_w_cuts_param.py --file /home/vault/capm/sn0527/PMT_LED_Temp/Run00705/fit_tree.root -q %f -s %f" %(9.0, 4.8)  )
#	os.system( "python /home/vault/capm/sn0527/PMT_LED_Temp/fit_1pe_w_cuts_param.py --file /home/vault/capm/sn0527/PMT_LED_Temp/Run00705/fit_tree.root -q %f -s %f" %(7.0, 4.9)  )
	
	
if __name__== "__main__":
	main()