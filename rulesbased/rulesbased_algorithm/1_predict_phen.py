# Predicts phenotype for combined RGI and EPI results 
# Creates a "isolatename.gp.txt" file for each isolate
# To use, put all RGI and EPI results (.txt) and phenotype.json into same directory as this file 
# python predict_phen.py phenotype.json


import sys
import csv
import json
import os
import glob

epi_csvs = glob.glob('*.effluxpumpmetamodels.txt')
rgi_csvs = glob.glob('*.output.txt')

for fr in rgi_csvs:
	rgi_dict = {}
	redundant = []
	sample= fr.split(".")[0]
	with open(fr, 'r') as rgi_file:
		rgi_file = csv.reader(rgi_file, delimiter='\t')
		header = next(rgi_file, None)
		for row in rgi_file:
			if row[21].isdigit():
				if row[5] == "Perfect" and row[21] not in rgi_dict.keys():
					rgi_dict.setdefault(row[21],[]).append(row[8])
					rgi_dict.setdefault(row[21],[]).append(row[5])
			elif row[22].isdigit():
				if row[5] == "Perfect" and row[22] not in rgi_dict.keys():
					rgi_dict.setdefault(row[22],[]).append(row[8])
					rgi_dict.setdefault(row[22],[]).append(row[5])				
	with open(fr, 'r') as rgi_file:
		rgi_file = csv.reader(rgi_file, delimiter='\t')
		header = next(rgi_file, None)
		for row in rgi_file:
			if row[21].isdigit():
				if row[5] == "Strict" and row[21] not in rgi_dict.keys():
					rgi_dict.setdefault(row[21],[]).append(row[8])
					rgi_dict.setdefault(row[21],[]).append(row[5])
			elif row[22].isdigit():
				if row[5] == "Strict" and row[22] not in rgi_dict.keys():
					rgi_dict.setdefault(row[22],[]).append(row[8])
					rgi_dict.setdefault(row[22],[]).append(row[5])
						
	for fe in epi_csvs:
		if sample in fe:
			perf_epi_metam = {}
			perf_epi_m = {}
			part_epi_metam = {}
			part_epi_m = {}
			put_epi_metam = {}
			put_epi_m = {}	
			with open(fe, 'r') as epi_file:
				epi_file = csv.reader(epi_file, delimiter='\t')
				header = next(epi_file, None)
				for row in epi_file:			
					if row[2] == "Perfect":
						perf_epi_metam.setdefault(row[0],[]).append(row[1])
						perf_epi_m.setdefault(row[3],[]).append(row[8])	
					elif row[2] == "Partial":
						part_epi_metam.setdefault(row[0],[]).append(row[1])
						part_epi_m.setdefault(row[3],[]).append(row[8])
					elif row[2] == "Putative":
						put_epi_metam.setdefault(row[0],[]).append(row[1])
						put_epi_m.setdefault(row[3],[]).append(row[8])	
	for fe in epi_csvs:
		if sample in fe:
			snp_epi_metam = {}
			with open(fe, 'r') as epi_file:
				epi_file = csv.reader(epi_file, delimiter='\t')
				header = next(epi_file, None)
				for row in epi_file:	
					if row[9] != "n/a":
						snp_epi_metam.setdefault(row[0],[]).append(row[1])

	for kr in rgi_dict.keys():		
		for ke in perf_epi_m.keys():
			if kr == ke:
				redundant.append(kr)
		for ke in part_epi_m.keys():
			if kr == ke:
				redundant.append(kr)		
		for ke in put_epi_m.keys():	
			if kr == ke:
				redundant.append(kr)	

	for x in set(redundant):
		del rgi_dict[x]

	with open(sys.argv[-1], 'r') as phen_file:
		with open(os.path.splitext(fr)[0]+".gp.txt", "w") as af:
			writer = csv.writer(af, delimiter='\t', dialect='excel')
			writer.writerow(["Meta-model ID / Model ID", "Meta-model name / Model name", "Category", "Confers resistance to"])																							
			phen_data = json.load(phen_file) 
			for entry in phen_data:
				phen = []
				for x in (phen_data[entry]["model_phenotypes"]):
					matchdict={}
					for k,v in snp_epi_metam.items():
						if k == entry:
							phen.append(phen_data[entry]["model_phenotypes"][x]["name"])
							matchdict[entry]=[k, v[0], "Partial", ", ".join(phen)]					
					for k,v in rgi_dict.items():
						if k == entry:
							phen.append(phen_data[entry]["model_phenotypes"][x]["name"])
							matchdict[entry]=[k, v[0], v[1], ", ".join(phen)]				
				for key, value in matchdict.items():
					writer.writerow(value)
