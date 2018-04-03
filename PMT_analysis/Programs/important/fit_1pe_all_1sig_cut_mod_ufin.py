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
	run_list = [ ["00695", "00696", "00697", "00698"], ["00700", "00701", "00702"], ["00703", "00704", "00705"], ["00706", "00707", "00708", "00709"] ]
	sigma1_limit_list = [ [ 2.9, 6.5, 6.5, 3.5 ], [ 3.4, 3.7, 6.5 ], [6.5, 3.7, 6.5], [ 3.4, 3.5, 6.5, 6.5 ] ]
	Q1_value = [ 7.0, 9.0, 9.0, 7.0 ]
	Temp_list = [ "+23", "-73", "-43", "-9"]

		
	# loop over temperatures
	#~ for j in range(4):
	for j in range(4):
		print "\n"
		print "--------------------------- new temperature ---------------------------"
		with open( "/home/vault/capm/sn0527/PMT_LED_Temp/Auswertung_ultrafinal/1sig_cut_mod/parameter_T=%sC_1sig_cut.txt" % (Temp_list[j]  ), "w" ) as plotfile:
			plotfile.write( "#entries\tchi2\tmu\tmu_error\tQ1\tQ1_error\tsigma_1\tsigma_1_error\tQ0\tQ0_error\tsigma_0\tsigma_0_error\tscaling\tscaling_error\tw\tw_error\talpha\talpha_error\n" )
			# loop over intensities
			for i in range(len(run_list[j])):
				print "\n"
				print "--------------------------- new Run ---------------------------"
				print "\n"
				parameters = ""
				os.system( "python /home/vault/capm/sn0527/PMT_LED_Temp/Auswertung_ultrafinal/1sig_cut_mod/fit_1pe_single_1sig_cut_mod_ufin.py --file /home/vault/capm/sn0527/PMT_LED_Temp/Run%s/fit_tree_1_sigma.root -q %f -s %f -c %s" %(  run_list[j][i], Q1_value[j], sigma1_limit_list[j][i], "1_sigma_cut" )  )
				os.system( "cp /home/vault/capm/sn0527/PMT_LED_Temp/Run%s/pmt_1pe_res_1sig_cut_unmod_final.pdf /home/vault/capm/sn0527/PMT_LED_Temp/Auswertung_ultrafinal/1sig_cut_mod/single_plots/pmt_1pe_1sig_lim_unmod_Run%s.pdf" %( run_list[j][i], run_list[j][i] ) )
				os.system( "cp /home/vault/capm/sn0527/PMT_LED_Temp/Run%s/fit_data.txt /home/vault/capm/sn0527/PMT_LED_Temp/Auswertung_ultrafinal/1sig_cut_mod/fit_data/fit_data_Run%s.txt" %( run_list[j][i], run_list[j][i] ) )
				os.system( "cp /home/vault/capm/sn0527/PMT_LED_Temp/Run%s/histogram_data.txt /home/vault/capm/sn0527/PMT_LED_Temp/Auswertung_ultrafinal/1sig_cut_mod/histogram_data/histogram_data_Run%s.txt" %( run_list[j][i], run_list[j][i] ) )
				os.system( "cp /home/vault/capm/sn0527/PMT_LED_Temp/Run%s/residual.txt /home/vault/capm/sn0527/PMT_LED_Temp/Auswertung_ultrafinal/1sig_cut_mod/residual/residual_Run%s.txt" %( run_list[j][i], run_list[j][i] ) )
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