import subprocess
import csv
import re
import matplotlib.pyplot as plt
import numpy as np
from downloadmatrix import *
import sys

def main():
	matrices = []
	#datafile = open('matrix_list_full.txt', 'r')
	datafile = open('matrix_list_revised.txt', 'r')
	myreader = csv.reader(datafile)
	hts = []
	mat_id = []
	DM = DownloadMatrix()
	libufgetpath = '/kuacc/users/nahmad16/.ufget/'

	for row in myreader:

		matrix_id = row[0]
		file_downloaded = 0
		while file_downloaded == 0:
			print('Downloading matrix {}'.format(matrix_id))
			[rc,path] = DM.downloadSSGET(matrix_id)
			if rc == 7:
			    file_downloaded = 0
			else:
			    file_downloaded = 1
		mtxfilenameindex = path.rindex("/")
		mtxfilenameindex = path[:mtxfilenameindex].rindex("/")

		ssgetfilepath = path[:path.rindex("/")] + ".tar.gz"
		libufgetsavepath = libufgetpath + path[36:mtxfilenameindex+1]
		result = subprocess.run(['mkdir', libufgetsavepath], stdout=subprocess.PIPE)
		result = subprocess.run(['cp', ssgetfilepath,libufgetsavepath], stdout=subprocess.PIPE)
		command = "./hts_test " + matrix_id + " build/ "
		process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
		output, error = process.communicate()
		str_output = output.decode('utf-8').split('\n')
		sptrsv_time = 0.0
		
		print(matrix_id)
		#for line in str_output:
		#print(str_output)
		for line in str_output:
			#print(line)	
			#text = re.findall(r'SpTRSV *[^\w ] *(.*)', line)
			#text1 = re.findall(r'Preprocess *[^\w ] *(.*)', line)
			text = re.findall(r'^SpTRSV *[^\w] *(.*)', line)
			text1 = re.findall(r'^Preprocess *[^\w] *(.*)', line)
			#print("Time" + str(text))
			#print("Preproces" + str(text1))
			
			#text = line.split(",")
			if len(text):
				print("Appending...")
				sptrsv_time = float(text[0])
			if len(text1):
				preprocess_time = float(text1[0])
				hts.append(matrix_id+","+str(sptrsv_time)+","+str(preprocess_time))
				#hts.append(sptrsv_time)
			
		
	#hybrid_vs_mkl = [float(i)/j for i,j in zip(mkl,hybrid)]
	#hybrid_vs_elmc = [float(i)/j for i,j in zip(elmc,hybrid)]
	#hybrid_speedup = np.minimum(hybrid_vs_elmc,hybrid_vs_mkl)

	with open('hts.txt','w') as hts_file:
		for item in hts:
			hts_file.write('%s\n'%item)


if __name__ == '__main__':
    main()
