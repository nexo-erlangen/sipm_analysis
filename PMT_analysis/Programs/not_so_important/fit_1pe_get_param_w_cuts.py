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
	#run_list=["00703", "00704", "00705"]
	run_list = [ ["00695", "00696", "00697", "00698"], ["00700", "00701", "00702"], ["00703", "00704", "00705"], ["00706", "00707", "00708", "00709"] ]
	#sigma1_limit_list = [5.5, 5.2, 4.8]
	sigma1_limit_list = [ [ 4.2, 4.3, 3.9, 4.6 ], [ 5.0, 5.3, 5.2 ], [5.0, 5.2, 4.8], [ 4.9, 4.8, 4.6, 4.5 ] ]
	Q1_value = [ 7.0, 9.0, 9.0, 7.0 ]
	Temp_list = [ "+23", "-73", "-43", "-9"]
	
	'''for j in range(4):
	#j = 3;
		result_list=[]
		for i in range(len( run_list[j])):
			print "Run:\t\t\t", run_list[j][i] 
			print "Q1 lower end:\t\t", Q1_value[j]
			print "sigma1 upper limit:\t", sigma1_limit_list[j][i]
			result_list.append(fit_histogram(run_list[j][i], Q1_value[j], sigma1_limit_list[j][i]))
		
		with open("/home/vault/capm/sn0527/PMT_LED_Temp/Auswertung_w_cuts/fit_parameters_T=" + Temp_list[j] + "C.dat", "w") as f:
			f.write("#mu\tmu_error\tQ1\tQ1_error\tsigma_1\tsigma_1_error\tQ0\tQ0_error\tsigma_0\tsigma_0_error\tscaling\tscaling_error\tw\tw_error\talpha\talpha_error\n")
			for result in result_list:
				for entry in result:
					f.write("%.5f\t" % entry)
				f.write("\n")
			f.write("\n")
	'''
	for j in range(4):
		for i in range(len(run_list[j])):
			os.system( "python /home/vault/capm/sn0527/PMT_LED_Temp/fit_1pe_w_cuts_param_res.py --file /home/vault/capm/sn0527/PMT_LED_Temp/Run%s/fit_tree.root -q %f -s %f" %(  run_list[j][i], Q1_value[j], sigma1_limit_list[j][i] )  )
			
#	os.system( "python /home/vault/capm/sn0527/PMT_LED_Temp/fit_1pe_w_cuts_param.py --file /home/vault/capm/sn0527/PMT_LED_Temp/Run00705/fit_tree.root -q %f -s %f" %(9.0, 4.8)  )
#	os.system( "python /home/vault/capm/sn0527/PMT_LED_Temp/fit_1pe_w_cuts_param.py --file /home/vault/capm/sn0527/PMT_LED_Temp/Run00705/fit_tree.root -q %f -s %f" %(7.0, 4.9)  )
	
	
if __name__== "__main__":
	main()