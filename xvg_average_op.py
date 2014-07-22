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
parser = argparse.ArgumentParser(prog = 'xvg_average_op', usage='', add_help = False, formatter_class = argparse.RawDescriptionHelpFormatter, description =\
'''
**********************************************
v''' + version_nb + '''
author: Jean Helie (jean.helie@bioch.ox.ac.uk)
git: https://github.com/jhelie/xvg_average_op
**********************************************

[ DESCRIPTION ]
 
This script calculate the average of order param data contained in several xvg files.

It also calculates the (unbiased) standard deviation by using the Bienaymé formula
to calculate the variance (http://en.wikipedia.org/wiki/Variance).

[ USAGE ]

Option	      Default  	Description                    
-----------------------------------------------------
-f			: xvg file(s)
-o		average	: name of outptut file
--membrane			: 'AM_zCter','AM_zNter','SMa','SMz' or 'POPC'
--comments	@,#	: lines starting with these characters will be considered as comment

Other options
-----------------------------------------------------
--version		: show version number and exit
-h, --help		: show this menu and exit
 
''')

#options
parser.add_argument('-f', nargs='+', dest='xvgfilenames', help=argparse.SUPPRESS, required=True)
parser.add_argument('--membrane', dest='membrane', choices=['AM_zCter','AM_zNter','SMa','SMz'], default='not specified', help=argparse.SUPPRESS, required=True)
parser.add_argument('--comments', nargs=1, dest='comments', default=['@,#'], help=argparse.SUPPRESS)

#other options
parser.add_argument('--version', action='version', version='%(prog)s v' + version_nb, help=argparse.SUPPRESS)
parser.add_argument('-h','--help', action='help', help=argparse.SUPPRESS)

#=========================================================================================
# store inputs
#=========================================================================================

args = parser.parse_args()
args.output_file = args.output_file[0]
args.nb_skipping = args.nb_skipping[0]
args.nb_smoothing = args.nb_smoothing[0]
args.nan2num = args.nan2num[0]

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

for f in args.xvgfilenames:
	if not os.path.isfile(f):
		print "Error: file " + str(f) + " not found."
		sys.exit(1)

if args.nan2num != "no":
	try:
		args.nan2num = float(args.nan2num)
	except:
		print "Error: --nan should be set to a float value."
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
	global first_col
	global files_columns
	global columns_names
	global label_xaxis
	global label_yaxis
	global weight_sum
	nb_rows = 0
	nb_cols = 0
	weight_sum = 0
	label_xaxis = "x axis"
	label_yaxis = "y axis"
	files_columns = {}
	columns_names = []
	
	for f_index in range(0,len(args.xvgfilenames)):
		progress = '\r -reading file ' + str(f_index+1) + '/' + str(len(args.xvgfilenames)) + '                      '  
		sys.stdout.flush()
		sys.stdout.write(progress)
		filename = args.xvgfilenames[f_index]
		tmp_nb_rows_to_skip = 0
		files_columns[filename] = {"leg2col": {}}
		files_columns[filename]["weight"] = 1
		#get file content
		with open(filename) as f:
			lines = f.readlines()
		
		#determine legends and nb of lines to skip
		for l_index in range(0,len(lines)):
			line = lines[l_index]
			if line[-1] == '\n':
				line = line[:-1]
			if line[0] in args.comments:
				tmp_nb_rows_to_skip += 1
				if "legend \"" in line:
					try:
						tmp_col = int(int(line.split("@ s")[1].split(" ")[0]) + 1)
						tmp_name = line.split("legend \"")[1][:-1]
						files_columns[filename]["leg2col"][tmp_name] = tmp_col
					except:
						print "\nError: unexpected data format in line " + str(l_index) + " in file " + str(filename) + "."
						print " -> " + str(line)
						sys.exit(1)
					if f_index == 0:
						if tmp_name in columns_names:
							print "znError: the legend '" + str(tmp_name) + "' is used twice in file " + str(filename) + "."
							sys.exit(1)
						else:
							columns_names.append(tmp_name)
					else:
						if tmp_name not in columns_names:
							print "\nError: legend '" + str(tmp_name) + "' is present in file " + str(filename) + " but not in " + str(args.xvgfilenames[0]) + "."
							sys.exit(1)
				if "xaxis" in line and  "label " in line:
					label_xaxis = line.split("label ")[1]
				if "yaxis" in line and  "label " in line:
					label_yaxis = line.split("label ")[1]
				if "weight" in line:
					if "-> weight = " in line:
						files_columns[filename]["weight"] = float(line.split("-> weight = ")[1])
						if files_columns[filename]["weight"] < 0:
							print "\nError: the weight in file " + str(filename) + " should be a positive number."
							print " -> " + str(line)
							sys.exit(1)
					else:
						print "\nWarning: keyword 'weight' found in the comments of file " + str(filename) + ", but weight not read in as the format '-> weight = ' wasn't found."
		
		#get data
		files_columns[filename]["data"] = np.loadtxt(filename, skiprows = tmp_nb_rows_to_skip)
		
		#check that each file has the same number of data rows
		if f_index == 0:
			nb_rows = np.shape(files_columns[filename]["data"])[0]
			if nb_rows < args.nb_smoothing:
				print "Error: the number of data rows (" + str(nb_rows) + ") is smaller than the option --smooth specified (" + str(args.nb_smoothing) + ")."
				sys.exit(1)
			if nb_rows < args.nb_skipping:
				print "Error: the number of data rows (" + str(nb_rows) + ") is smaller than the option --skip specified (" + str(args.nb_skipping) + ")."
				sys.exit(1)
		else:
			if np.shape(files_columns[filename]["data"])[0] != nb_rows:
				print "Error: file " + str(filename) + " has " + str(np.shape(files_columns[filename]["data"])[0]) + " data rows, whereas file " + str(args.xvgfilenames[0]) + " has " + str(nb_rows) + " data rows."
				sys.exit(1)
	
		#check that each file has the same number of columns
		if f_index == 0:
			nb_cols = np.shape(files_columns[filename]["data"])[1]
		else:
			if np.shape(files_columns[filename]["data"])[1] != nb_cols:
				print "Error: file " + str(filename) + " has " + str(np.shape(files_columns[filename]["data"])[1]) + " data columns, whereas file " + str(args.xvgfilenames[0]) + " has " + str(nb_cols) + " data columns."
				sys.exit(1)
			
		#check that each file has the same first column
		if f_index == 0:
			first_col = files_columns[filename]["data"][:,0]
		else:
			if not np.array_equal(files_columns[filename]["data"][:,0],first_col):
				print "\nError: the first column of file " + str(filename) + " is different than that of " + str(args.xvgfilenames[0]) + "."
				sys.exit(1)
		
		#update weight sum
		weight_sum += files_columns[filename]["weight"]

	return

#=========================================================================================
# core functions
#=========================================================================================

def rolling_avg(loc_list):												#DONE

	loc_arr = np.asarray(loc_list)
	shape = (loc_arr.shape[-1]-args.nb_smoothing+1,args.nb_smoothing)
	strides = (loc_arr.strides[-1],loc_arr.strides[-1])   	
	return scipy.stats.nanmean(np.lib.stride_tricks.as_strided(loc_arr, shape=shape, strides=strides), -1)
def calculate_avg():													#DONE

	global data_avg
	global data_std
	global nb_rows
	global nb_cols
	
	#calculate raw average
	#---------------------
	data_avg = np.zeros((nb_rows,nb_cols))
	if len(args.xvgfilenames) > 1:
		data_std = np.zeros((nb_rows,nb_cols-1))
	data_avg[:,0] = first_col
	for col_index in range(1, nb_cols):
		col_name = columns_names[col_index-1]
		#initialise average with first file
		filename = args.xvgfilenames[0]
		tmp_col_nb = files_columns[filename]["leg2col"][col_name]
		tmp_col_avg = files_columns[filename]["data"][:,tmp_col_nb:tmp_col_nb+1] * files_columns[filename]["weight"] * len(args.xvgfilenames) / float(weight_sum)
		
		#add columns of following files
		for f_index in range(1,len(args.xvgfilenames)):
			filename = args.xvgfilenames[f_index]
			tmp_col_nb = files_columns[filename]["leg2col"][col_name]
			tmp_col_avg = np.concatenate([tmp_col_avg,files_columns[filename]["data"][:,tmp_col_nb:tmp_col_nb+1] * files_columns[filename]["weight"] * len(args.xvgfilenames) / float(weight_sum)], axis = 1)	
				
		if len(args.xvgfilenames) > 1:

			#calculate weighted average taking into account "nan"
			#----------------------------------------------------
			data_avg[:,col_index] =  scipy.stats.nanmean(tmp_col_avg, axis = 1)
						
			#calculate unbiased weighted std dev taking into account "nan"
			#-------------------------------------------------------------
			#initialise average with first file
			filename = args.xvgfilenames[0]
			tmp_col_nb = files_columns[filename]["leg2col"][col_name]
			tmp_col_std = files_columns[filename]["weight"] * (files_columns[filename]["data"][:,tmp_col_nb:tmp_col_nb+1] - data_avg[:,col_index:col_index+1])**2
			tmp_weigh2_sum = files_columns[filename]["weight"]**2
			
			#add columns of following files
			for f_index in range(1,len(args.xvgfilenames)):
				filename = args.xvgfilenames[f_index]
				tmp_col_nb = files_columns[filename]["leg2col"][col_name]
				tmp_col_std = np.concatenate([tmp_col_std, files_columns[filename]["weight"] * (files_columns[filename]["data"][:,tmp_col_nb:tmp_col_nb+1] - data_avg[:,col_index:col_index+1])**2], axis = 1)	
				tmp_weigh2_sum += files_columns[filename]["weight"]**2
						
			#calculate unbiased standard deviation as defined on wikipedia: https://en.wikipedia.org/wiki/Weighted_variance#Weighted_sample_variance
			tmp_col_std = np.sqrt(weight_sum / float(weight_sum**2 - tmp_weigh2_sum) * scipy.nansum(tmp_col_std, axis = 1))
			data_std[:,col_index-1] = tmp_col_std

		else:
			data_avg[:,col_index] = tmp_col_avg[:,0]
			
	#update by smoothing
	#-------------------
	if args.nb_smoothing > 1:
		nb_rows = nb_rows - args.nb_smoothing + 1
		tmp_data_avg_smoothed = np.zeros((nb_rows,nb_cols))
		if len(args.xvgfilenames) > 1:
			tmp_data_std_smoothed = np.zeros((nb_rows,nb_cols-1))
		tmp_data_avg_smoothed[:,0] = np.transpose(rolling_avg(np.transpose(data_avg[:,0])))
		for col_index in range(1, nb_cols):
			tmp_data_avg_smoothed[:,col_index] = np.transpose(rolling_avg(np.transpose(data_avg[:,col_index])))
			if len(args.xvgfilenames) > 1:
				tmp_data_std_smoothed[:,col_index-1] = np.transpose(rolling_avg(np.transpose(data_std[:,col_index-1])))
		data_avg = tmp_data_avg_smoothed
		if len(args.xvgfilenames) > 1:
			data_std = tmp_data_std_smoothed
	
	#update by skipping
	#------------------
	if args.nb_skipping > 1 :
		rows_to_keep = [r for r in range(0,nb_rows) if r%args.nb_skipping ==0]
		nb_rows = len(rows_to_keep)
		data_avg = data_avg[rows_to_keep,:]
		if len(args.xvgfilenames) > 1:
			data_std = data_std[rows_to_keep,:]
	
	#replace nan values if necessary
	#-------------------------------
	if args.nan2num != "no":
		data_avg[np.isnan(data_avg)] = args.nan2num
		if len(args.xvgfilenames) > 1:
			data_std[np.isnan(data_std)] = args.nan2num
	
	return

#=========================================================================================
# outputs
#=========================================================================================

def write_xvg():														#DONE

	#open files
	filename_xvg = os.getcwd() + '/' + str(args.output_file) + '.xvg'
	output_xvg = open(filename_xvg, 'w')
	
	#general header
	output_xvg.write("# [average xvg - written by xvg_average v" + str(version_nb) + "]\n")
	tmp_files = ""
	for f in args.xvgfilenames:
		tmp_files += "," + str(f)
	output_xvg.write("# - files: " + str(tmp_files[1:]) + "\n")
	output_xvg.write("# - skipping: " + str(args.nb_skipping) + " frames\n")
	output_xvg.write("# - smoothing: " + str(args.nb_smoothing) + " frames\n")
	if weight_sum > len(args.xvgfilenames):
		output_xvg.write("# -> weight = " + str(weight_sum) + "\n")
	
	#xvg metadata
	output_xvg.write("@ title \"Average xvg\"\n")
	output_xvg.write("@ xaxis label " + str(label_xaxis) + "\n")
	output_xvg.write("@ yaxis label " + str(label_yaxis) + "\n")
	output_xvg.write("@ autoscale ONREAD xaxes\n")
	output_xvg.write("@ TYPE XY\n")
	output_xvg.write("@ view 0.15, 0.15, 0.95, 0.85\n")
	output_xvg.write("@ legend on\n")
	output_xvg.write("@ legend box on\n")
	output_xvg.write("@ legend loctype view\n")
	output_xvg.write("@ legend 0.98, 0.8\n")
	if len(args.xvgfilenames) > 1:
		output_xvg.write("@ legend length " + str((nb_cols-1)*2) + "\n")
		for col_index in range(0,nb_cols-1):
			output_xvg.write("@ s" + str(col_index) + " legend \"" + str(columns_names[col_index]) + " (avg)\"\n")
		for col_index in range(0,nb_cols-1):
			output_xvg.write("@ s" + str(nb_cols - 1 + col_index) + " legend \"" + str(columns_names[col_index]) + " (std)\"\n")
	else:
		output_xvg.write("@ legend length " + str(nb_cols-1) + "\n")
		for col_index in range(0,nb_cols-1):
			output_xvg.write("@ s" + str(col_index) + " legend \"" + str(columns_names[col_index]) + "\"\n")
	
	#data
	for r in range(0, nb_rows):
		results = str(data_avg[r,0])
		#avg
		for col_index in range(1,nb_cols):
			results += "	" + "{:.6e}".format(data_avg[r,col_index])
		#std
		if len(args.xvgfilenames) > 1:
			for col_index in range(0,nb_cols-1):
				results += "	" + "{:.6e}".format(data_std[r,col_index])
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
