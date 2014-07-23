#generic python modules
import argparse
import operator
from operator import itemgetter
import sys, os, shutil
import os.path

##########################################################################################
# RETRIEVE USER INPUTS
##########################################################################################

#=========================================================================================
# create parser
#=========================================================================================
version_nb = "0.0.1"
parser = argparse.ArgumentParser(prog = 'xvg_average_op_complex', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/xvg_average_op
**********************************************

[ DESCRIPTION ]
 
This script calculate the average of order param data contained in several xvg files.

It also calculates the (unbiased) standard deviation by using the Bienayme formula
to calculate the variance (http://en.wikipedia.org/wiki/Variance).

NB:
the script may give out a warning 'return np.mean(x,axis)/factor', it's ok. it's just
scipy warning us that there were only nans on a row, the result will be a nan as we
expect (see this thread: https://github.com/scipy/scipy/issues/2898).

[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: xvg file(s)
-o		op_avg	: name of outptut file
--membrane		: 'AM_zCter','AM_zNter','SMa','SMz' or 'POPC'
--comments	@,#	: lines starting with these characters will be considered as comment

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#options
parser.add_argument('-f', nargs='+', dest='xvgfilenames', help=argparse.SUPPRESS, required=True)
parser.add_argument('-o', nargs=1, dest='output_file', default=["op_avg"], help=argparse.SUPPRESS)
parser.add_argument('--membrane', dest='membrane', choices=['AM_zCter','AM_zNter','SMa','SMz','POPC'], default='not specified', help=argparse.SUPPRESS, required=True)
parser.add_argument('--comments', nargs=1, dest='comments', default=['@,#'], help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
args.output_file = args.output_file[0]
args.comments = args.comments[0].split(',')

#=========================================================================================
# import modules (doing it now otherwise might crash before we can display the help menu!)
#=========================================================================================

#generic science modules
try:
	import numpy as np
except:
	print "Error: you need to install the np module."
	sys.exit(1)
try:
	import scipy
	import scipy.stats
except:
	print "Error: you need to install the scipy module."
	sys.exit(1)

#=======================================================================
# sanity check
#=======================================================================

if len(args.xvgfilenames) == 1:
	print "Error: only 1 data file specified."
	sys.exit(1)
	
for f in args.xvgfilenames:
	if not os.path.isfile(f):
		print "Error: file " + str(f) + " not found."
		sys.exit(1)

##########################################################################################
# FUNCTIONS DEFINITIONS
##########################################################################################

#=========================================================================================
# data loading
#=========================================================================================

def load_xvg():															#DONE
	
	global nb_rows
	global nb_cols
	global weights
	global data_op_upper_avg
	global data_op_upper_std
	global data_op_upper_nb
	global data_op_lower_avg
	global data_op_lower_std
	global data_op_lower_nb
	nb_rows = 0
	nb_cols = 0
	weights = np.ones(len(args.xvgfilenames))
		
	for f_index in range(0,len(args.xvgfilenames)):
		#display progress
		progress = '\r -reading file ' + str(f_index+1) + '/' + str(len(args.xvgfilenames)) + '                      '  
		sys.stdout.flush()
		sys.stdout.write(progress)
		
		#get file content
		filename = args.xvgfilenames[f_index]
		with open(filename) as f:
			lines = f.readlines()
		
		#determine legends and nb of lines to skip
		tmp_nb_rows_to_skip = 0
		for l_index in range(0,len(lines)):
			line = lines[l_index]
			if line[0] in args.comments:
				tmp_nb_rows_to_skip += 1
				if "weight" in line:
					if "-> weight = " in line:
						weights[f_index] = float(line.split("-> weight = ")[1])
						if weights[f_index] < 0:
							print "\nError: the weight in file " + str(filename) + " should be a positive number."
							print " -> " + str(line)
							sys.exit(1)
					else:
						print "\nWarning: keyword 'weight' found in the comments of file " + str(filename) + ", but weight not read in as the format '-> weight = ' wasn't found."
		
		#get data
		tmp_data = np.loadtxt(filename, skiprows = tmp_nb_rows_to_skip)
		
		#check that each file has the same number of data rows
		if f_index == 0:
			nb_rows = np.shape(tmp_data)[0]
			data_op_upper_avg = np.zeros((nb_rows, len(args.xvgfilenames) + 1))			#distance, avg op upper for each file
			data_op_upper_std = np.zeros((nb_rows, len(args.xvgfilenames)))				#std op upper for each file
			data_op_upper_nb = np.zeros((nb_rows, len(args.xvgfilenames)))				#nb op upper for each file
			data_op_lower_avg = np.zeros((nb_rows, len(args.xvgfilenames) + 1))			#distance, avg op upper for each file
			data_op_lower_std = np.zeros((nb_rows, len(args.xvgfilenames)))				#std op upper for each file
			data_op_lower_nb = np.zeros((nb_rows, len(args.xvgfilenames)))				#nb op upper for each file
		else:
			if np.shape(tmp_data)[0] != nb_rows:
				print "Error: file " + str(filename) + " has " + str(np.shape(tmp_data)[0]) + " data rows, whereas file " + str(args.xvgfilenames[0]) + " has " + str(nb_rows) + " data rows."
				sys.exit(1)
		#check that each file has the same number of columns
		if f_index == 0:
			nb_cols = np.shape(tmp_data)[1]
		else:
			if np.shape(tmp_data)[1] != nb_cols:
				print "Error: file " + str(filename) + " has " + str(np.shape(tmp_data)[1]) + " data columns, whereas file " + str(args.xvgfilenames[0]) + " has " + str(nb_cols) + " data columns."
				sys.exit(1)
		#check that each file has the same first column
		if f_index == 0:
			data_op_upper_avg[:,0] = tmp_data[:,0]
			data_op_lower_avg[:,0] = tmp_data[:,0]
		else:
			if not np.array_equal(tmp_data[:,0],data_op_upper_avg[:,0]):
				print "\nError: the first column of file " + str(filename) + " is different than that of " + str(args.xvgfilenames[0]) + "."
				sys.exit(1)
		
		#store data
		if args.membrane == "AM_zCter":
			data_op_upper_avg[:, f_index + 1] = tmp_data[:,3]
			data_op_upper_std[:, f_index] = tmp_data[:,6]
			data_op_upper_nb[:, f_index] = tmp_data[:,9]
			data_op_lower_avg[:, f_index + 1] = tmp_data[:,13]
			data_op_lower_std[:, f_index] = tmp_data[:,17]
			data_op_lower_nb[:, f_index] = tmp_data[:,21]
		elif args.membrane == "AM_zNter":
			data_op_upper_avg[:, f_index + 1] = tmp_data[:,15]
			data_op_upper_std[:, f_index] = tmp_data[:,18]
			data_op_upper_nb[:, f_index] = tmp_data[:,21]
			data_op_lower_avg[:, f_index + 1] = tmp_data[:,4]
			data_op_lower_std[:, f_index] = tmp_data[:,8]
			data_op_lower_nb[:, f_index] = tmp_data[:,12]
		elif args.membrane == "SMa":
			data_op_upper_avg[:, f_index + 1] = tmp_data[:,16]
			data_op_upper_std[:, f_index] = tmp_data[:,20]
			data_op_upper_nb[:, f_index] = tmp_data[:,24]
			data_op_lower_avg[:, f_index + 1] = tmp_data[:,4]
			data_op_lower_std[:, f_index] = tmp_data[:,8]
			data_op_lower_nb[:, f_index] = tmp_data[:,12]
		elif args.membrane == "SMz":
			data_op_upper_avg[:, f_index + 1] = tmp_data[:,12]
			data_op_upper_std[:, f_index] = tmp_data[:,15]
			data_op_upper_nb[:, f_index] = tmp_data[:,18]
			data_op_lower_avg[:, f_index + 1] = tmp_data[:,3]
			data_op_lower_std[:, f_index] = tmp_data[:,6]
			data_op_lower_nb[:, f_index] = tmp_data[:,9]
		elif args.membrane == "POPC":
			data_op_upper_avg[:, f_index + 1] = tmp_data[:,8]
			data_op_upper_std[:, f_index] = tmp_data[:,8]
			data_op_upper_nb[:, f_index] = tmp_data[:,12]
			data_op_lower_avg[:, f_index + 1] = tmp_data[:,2]
			data_op_lower_std[:, f_index] = tmp_data[:,4]
			data_op_lower_nb[:, f_index] = tmp_data[:,6]
	return

#=========================================================================================
# core functions
#=========================================================================================

def calculate_avg():													#DONE

	global avg_op_upper_avg
	global avg_op_upper_std
	global avg_op_lower_avg
	global avg_op_lower_std
				
	avg_op_upper_avg = np.zeros((nb_rows, 2))
	avg_op_upper_std = np.zeros((nb_rows, 1))
	avg_op_lower_avg = np.zeros((nb_rows, 2))
	avg_op_lower_std = np.zeros((nb_rows, 1))

	#distances
	avg_op_upper_avg[:,0] = data_op_upper_avg[:,0]
	avg_op_lower_avg[:,0] = data_op_lower_avg[:,0]

	#calculate weighted average taking into account "nan"
	#----------------------------------------------------
	avg_op_upper_avg[:,1] =  scipy.stats.nanmean(data_op_upper_avg[:,1:] * weights * len(args.xvgfilenames) / float(np.sum(weights)) , axis = 1)
	avg_op_lower_avg[:,1] =  scipy.stats.nanmean(data_op_lower_avg[:,1:] * weights * len(args.xvgfilenames) / float(np.sum(weights)) , axis = 1)

	#calculate unbiased weighted std dev taking into account "nan"
	#-------------------------------------------------------------
	#from Bienayme formula
	# var(Xavg) = 1/(sum(wi))**2 * sum(wi**2 * var(Xi))
		
	#calculate total number of points
	tmp_nb_total_upper = np.copy(data_op_upper_nb)
	tmp_nb_total_upper[tmp_nb_total_upper != 0] += 1
	tmp_nb_total_upper = np.sum(tmp_nb_total_upper, axis = 1)
	tmp_nb_total_upper -= 1
	tmp_nb_total_upper[tmp_nb_total_upper == 0] = 1
	tmp_nb_total_upper[tmp_nb_total_upper == -1] = 1
	
	#calculate total number of points
	tmp_nb_total_lower = np.copy(data_op_lower_nb)
	tmp_nb_total_lower[tmp_nb_total_lower != 0] += 1
	tmp_nb_total_lower = np.sum(tmp_nb_total_lower, axis = 1)
	tmp_nb_total_lower -= 1
	tmp_nb_total_lower[tmp_nb_total_lower == 0] = 1
	tmp_nb_total_lower[tmp_nb_total_lower == -1] = 1
	
	#apply bienayme formula
	avg_op_upper_std[:,0] = np.sqrt(np.nansum(weights**2 * data_op_upper_std**2 * data_op_upper_nb, axis = 1) / (np.sum(weights)**2 * tmp_nb_total_upper))
	avg_op_lower_std[:,0] = np.sqrt(np.nansum(weights**2 * data_op_lower_std**2 * data_op_lower_nb, axis = 1) / (np.sum(weights)**2 * tmp_nb_total_lower))
		
	return

#=========================================================================================
# outputs
#=========================================================================================

def write_xvg():														#DONE

	#open files
	filename_xvg = os.getcwd() + '/' + str(args.output_file) + '.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [average xvg - written by xvg_average_op v" + str(version_nb) + "]\n")
	tmp_files = ""
	for f in args.xvgfilenames:
		tmp_files += "," + str(f)
	output_xvg.write("# - files: " + str(tmp_files[1:]) + "\n")
	if np.sum(weights) > len(args.xvgfilenames):
		output_xvg.write("# -> weight = " + str(np.sum(weights)) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Average xvg\"\n")
	output_xvg.write("@ xaxis label \"distance from cluster z axis (Angstrom)\"\n")
	output_xvg.write("@ yaxis label \"order parameter\"\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	output_xvg.write("@ legend length 4\n")
	output_xvg.write("@ s0 legend \"upper (avg)\"\n")
	output_xvg.write("@ s1 legend \"upper (std)\"\n")
	output_xvg.write("@ s2 legend \"lower (avg)\"\n")
	output_xvg.write("@ s3 legend \"lower (std)\"\n")
	
	#data
	for r in range(0, nb_rows):
		results = str(avg_op_upper_avg[r,0])
		results += "	" + "{:.6e}".format(avg_op_upper_avg[r,1]) + "	" + "{:.6e}".format(avg_op_upper_std[r,0]) + "	" + "{:.6e}".format(avg_op_lower_avg[r,1]) + "	" + "{:.6e}".format(avg_op_lower_std[r,0])
		output_xvg.write(results + "\n")		
	output_xvg.close()	
	
	return

##########################################################################################
# MAIN
##########################################################################################

print "\nReading files..."
load_xvg()

print "\n\nWriting average file..."
calculate_avg()
write_xvg()

#=========================================================================================
# exit
#=========================================================================================
print "\nFinished successfully! Check result in file '" + args.output_file + ".xvg'."
print ""
sys.exit(0)
