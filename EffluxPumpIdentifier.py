import os
import sys
import logging
import json

class EffluxPump(object):
	"""Class for efflux pump searches."""
	def __init__(self, filepath):
		self.filepath = filepath

	def __repr__(self):
		"""Returns EffluxPump class full object."""
		return "EffluxPump({}".format(self.__dict__)

	def run(self):

		"""Parses CARD.json file for efflux pump meta model info."""
		metamodelname = {}
		metamodelid = {}
		components = {}
		desc = {}
		effluxpumpmetamodels ={}
		
		with open("/workspace/tsangkk2/card.json") as card_file:
			card_data = json.load(card_file)

		del card_data["_version"]
		del card_data["_comment"]
		del card_data["_timestamp"]

		for entry in card_data: 
			if card_data[entry]["model_type_id"] == "41112":
				metamodelid[card_data[entry]["ARO_name"]] = card_data[entry]["ARO_id"]
				desc[card_data[entry]["ARO_name"]] = card_data[entry]["ARO_description"]
				try:
					components[entry] = card_data[entry]["model_param"]["41141"]["param_value"][list(card_data[entry]["model_param"]["41141"]["param_value"])[0]].split(',')
				except:
					pass
			else:
				continue

		"""Parses RGI results JSON file for efflux pump meta model results."""

		regulator = {}
		subunit = {}
		s_rgicardmatch = {}
		s_metamatch = {}
		r_rgicardmatch = {}
		r_metamatch = {}
		perfect_r = {}
		perfect_s = {}
		perfect_s1 = {}
		perfect_r1 = {}
		strict_r = {}
		strict_s = {}
		loose_r = {}
		loose_s ={}
		strict_r1 = {}
		strict_s1 = {}
		loose_r1 = {}
		loose_s1 ={}
		meta_match = {}
		perfect_metamatch = {}
		perfect_metamatch2 = {}
		partial_metamatch = {}
		partial_metamatch2 = {}
		partial_metamatch3 = {}
		partial_metamatch4 = {}

		with open (self.filepath) as rgi_file:
			rgi_data = json.load(rgi_file)

		try:
			del rgi_data["_metadata"]
		except:
			pass

		oldrgi_filename=str(rgi_file)
		newrgi_filename=oldrgi_filename.replace("<_io.TextIOWrapper name='","").replace("' mode='r' encoding='UTF-8'>","")
#		logger.info("WORKING ON "+newrgi_filename)

		"""Creates dictionaries for RGI hits that are efflux regulators and subunits."""

		for hsp in rgi_data:
			if (len(rgi_data[hsp])) >= 1:
				for hit in rgi_data[hsp]:
					if "36590" in rgi_data[hsp][hit]["ARO_category"]:
						regulator[hit] = rgi_data[hsp][hit]
						for metamodel,modellist in components.items():
							if regulator[hit]["model_id"] in modellist:
								x = regulator[hit]["model_id"]
								r_rgicardmatch.setdefault(x,[]).append(metamodel)
						for k,v in r_rgicardmatch.items():
							for model in v:
								for x in r_rgicardmatch.values():
									if model in x:
										r_metamatch.setdefault(model,[]).append(k)
										r_metamatch = dict((k, list(set(v))) for k, v in r_metamatch.items())
					elif "36298" in rgi_data[hsp][hit]["ARO_category"]:
						subunit[hit] = rgi_data[hsp][hit]
						for metamodel,modellist in components.items():
							if subunit[hit]["model_id"] in modellist:
								x = subunit[hit]["model_id"]
								s_rgicardmatch.setdefault(x,[]).append(metamodel)
						for k,v in s_rgicardmatch.items():
							for model in v:
								for x in s_rgicardmatch.values():
									if model in x:
										s_metamatch.setdefault(model,[]).append(k)
										s_metamatch = dict((k, list(set(v))) for k, v in s_metamatch.items())

			"""Creates dictionaries for perfect, strict, and loose subunits and regulators with query name as key."""
			for k,v in r_metamatch.items():
				for model in v:
					if (len(rgi_data[hsp])) >= 1:
						for x in rgi_data[hsp].values():						
							if model == x["model_id"] and x["type_match"] == "Perfect":
								perfect_r.setdefault(hsp, []).append((k, x))
							elif model == x["model_id"] and x["type_match"] == "Strict":
								if "snp" in rgi_data[hsp][hit]:
									perfect_r.setdefault(hsp, []).append((k, x))
								else:
									strict_r.setdefault(hsp, []).append((k, x))
							elif model == x["model_id"] and x["type_match"] == "Loose":
								loose_r.setdefault(hsp, []).append((k, x))
			for k,v in s_metamatch.items():
				for model in v:
					if (len(rgi_data[hsp])) >= 1:
						for x in rgi_data[hsp].values():
							if model == x["model_id"] and x["type_match"] == "Perfect":
								perfect_s.setdefault(hsp, []).append((k, x))
							elif model == x["model_id"] and x["type_match"] == "Strict":
								if "snp" in rgi_data[hsp][hit]:
									perfect_s.setdefault(hsp, []).append((k, x))
								else:
									strict_s.setdefault(hsp, []).append((k, x))
							elif model == x["model_id"] and x["type_match"] == "Loose":
								loose_s.setdefault(hsp, []).append((k, x))

		"""Creates dictionaries for perfect subunits and regulators with meta-model ID as the key."""
		for value in perfect_r.values():
			for x in value:
				for key in components.keys():
					if x[0] is key:
						perfect_r1.setdefault(key, []).append(x)

		for value in perfect_s.values():
			for x in value:
				for key in components.keys():
					if x[0] is key:
						perfect_s1.setdefault(key, []).append(x)

		"""Creates dictionaries for strict subunits and regulators with meta-model ID as the key."""
		for value in strict_r.values():
			for x in value:
				for key in components.keys():
					if x[0] is key:
						strict_r1.setdefault(key, []).append(x)

		for value in strict_s.values():
			for x in value:
				for key in components.keys():
					if x[0] is key:
						strict_s1.setdefault(key, []).append(x)

		"""Creates dictionaries for loose subunits and regulators with meta-model ID as the key."""
		for value in loose_r.values():
			for x in value:
				for key in components.keys():
					if x[0] is key:
						loose_r1.setdefault(key, []).append(x)

		for value in loose_s.values():
			for x in value:
				for key in components.keys():
					if x[0] is key:
						loose_s1.setdefault(key, []).append(x)

		"""Sort all perfect, strict, and loose hits with corresponding meta-models (key)"""
		
		for k1,v1 in perfect_s1.items():
			if type(v1) is list:
				for x in v1:
					meta_match.setdefault(k1, []).append(x)
			else:
				meta_match.setdefault(k1, []).append(v1)
		
		for k1,v1 in perfect_r1.items():
			if type(v1) is list:
				for x in v1:
					if x not in meta_match.values():
						meta_match.setdefault(k1, []).append(x)
			else:
				meta_match.setdefault(k1, []).append(v1)

		for k1,v1 in strict_s1.items():
			if type(v1) is list:
				for x in v1:
					if x not in meta_match.values():
						meta_match.setdefault(k1, []).append(x)
			else:
				meta_match.setdefault(k1, []).append(v1)

		for k1,v1 in strict_r1.items():		
			if type(v1) is list:
				for x in v1:
					if x not in meta_match.values():
						meta_match.setdefault(k1, []).append(x)
			else:
				meta_match.setdefault(k1, []).append(v1)

		for k1,v1 in loose_s1.items():
			if type(v1) is list:
				for x in v1:
					if x not in meta_match.values():
						meta_match.setdefault(k1, []).append(x)
			else:
				meta_match.setdefault(k1, []).append(v1)

		for k1,v1 in loose_r1.items():
			if type(v1) is list:
				for x in v1:
					if x not in meta_match.values():
						meta_match.setdefault(k1, []).append(x)
			else:
				meta_match.setdefault(k1, []).append(v1)

		"""PERFECT Dictionary with perfect and snp hits"""
		for k, v in meta_match.items():
			for k1,v1 in components.items():
				if k == k1:
					for metamodel_hit in v:
						for model in metamodel_hit:
							if type(model) is dict:
								for y in v1: #each model in meta-model in CARD
									if model["model_id"] == y and model["type_match"] == "Perfect":
										perfect_metamatch.setdefault(k, []).append(model)
									if model["model_id"] == y and model["type_match"] == "Strict" and "snp" in model:
										perfect_metamatch.setdefault(k, []).append(model)

		"""Creates PERFECT dictionary without duplicate values"""
		for k, v in perfect_metamatch.items():
			for i in range(0, len(v)):
				if v[i] not in v[i+1:]:
					perfect_metamatch2.setdefault(k, []).append(v[i])

		"""PERFECT efflux pump meta-models with snp"""	

		for key, value in perfect_metamatch2.items():
			perfect_jsondict ={}
			x3 = {}
			x2 = {}
			bitscore_dict ={}
			percidentity_dict ={}			
			perf_modellist =[]
			for k1,v1 in components.items():
				if key == k1:
					for x in value:
						newdict = {}
						perf_modellist.append(x["model_id"])
						perf_modellist1 = set(perf_modellist)
						if all(y in perf_modellist1 for y in v1):
							if any("snp" in x for x in value):			
								newdict = perfect_metamatch2[key]
								count=0
								for x in newdict:
									if x["type_match"] == "Perfect" and ("model_id: "+x["model_id"]) not in x2.keys():
										bitscore_dict.setdefault((x["model_id"]), []).append(x["bit_score"])
										bitscore_dict = dict((k, list(set(v))) for k, v in bitscore_dict.items())
										percidentity_dict.setdefault((x["model_id"]), []).append(x["perc_identity"])
										percidentity_dict = dict((k, list(set(v))) for k, v in percidentity_dict.items())
										for k,v in bitscore_dict.items():
											if x["bit_score"] == max(v):
												for k,v in percidentity_dict.items():
													if x["perc_identity"] == max(v):
														x2 = {"model_id: "+x["model_id"]: x}
									elif x["type_match"] == "Strict" and ("model_id: "+x["model_id"]) not in x2.keys():
										if "snp" in x:				
											if x not in x2.values():
												count+=1
												x2 = {"model_id: "+x["model_id"]+"_"+str(count): x}	
									x3.update(x2)
								perfect_jsondict["Perfect"] = x3
								for entry in card_data:
									if card_data[entry]["model_id"] == key:							
										effluxpumpmetamodels["meta-model ID: "+key+", "+card_data[entry]["model_name"]] = perfect_jsondict

		"""PERFECT efflux pump meta-models with no snp"""	
		for key, value in perfect_metamatch2.items():
			perfect_jsondict ={}
			x3 = {}
			x2 = {}
			bitscore_dict ={}
			percidentity_dict ={}			
			perf_modellist =[]
			for k1,v1 in components.items():
				if key == k1:
					for x in value:
						newdict = {}
						perf_modellist.append(x["model_id"])
						perf_modellist1 = set(perf_modellist)
						if all(y in perf_modellist1 for y in v1):
							if all("snp" not in x for x in value):			
								newdict = perfect_metamatch2[key]
								for x in newdict:
									if x["type_match"] == "Perfect" and ("model_id: "+x["model_id"]) not in x2.keys():
										bitscore_dict.setdefault((x["model_id"]), []).append(x["bit_score"])
										bitscore_dict = dict((k, list(set(v))) for k, v in bitscore_dict.items())
										percidentity_dict.setdefault((x["model_id"]), []).append(x["perc_identity"])
										percidentity_dict = dict((k, list(set(v))) for k, v in percidentity_dict.items())
										for k,v in bitscore_dict.items():
											if x["bit_score"] == max(v):
												for k,v in percidentity_dict.items():
													if x["perc_identity"] == max(v):
														x2 = {"model_id: "+x["model_id"]: x}
									x3.update(x2)
								perfect_jsondict["Perfect"] = x3
								for entry in card_data:
									if card_data[entry]["model_id"] == key:							
										effluxpumpmetamodels["meta-model ID: "+key+", "+card_data[entry]["model_name"]] = perfect_jsondict

		"""PARTIAL Dictionaries, one with perfect & snp hits, other with strict and loose hits"""
		for k, v in meta_match.items():
			partialmodel =[]
			for k1,v1 in components.items():
				if k == k1:
					for metamodel_hit in v:
						for model in metamodel_hit:
							if type(model) is dict:
								for y in v1: #each model in meta-model in CARD
									if model["model_id"] == y and model['type_match'] == "Perfect":
										partial_metamatch.setdefault(k, []).append(model)
									elif model["model_id"] == y and model['type_match'] == "Strict" and "snp" in x:
										partial_metamatch.setdefault(k, []).append(model)
									if partial_metamatch:
										for key, value in partial_metamatch.items():
											for x in value:
												partialmodel.append(x["model_id"])
											if model["model_id"] == y and model['type_match'] == "Strict" and model["model_id"] not in partialmodel:
												partial_metamatch2.setdefault(k, []).append(model)
											if model["model_id"] == y and model['type_match'] == "Loose" and model["model_id"] not in partialmodel:
												partial_metamatch2.setdefault(k, []).append(model)
									elif not partial_metamatch:
										if model["model_id"] == y and model['type_match'] == "Strict":
											partial_metamatch2.setdefault(k, []).append(model)
										if model["model_id"] == y and model['type_match'] == "Loose":
											partial_metamatch2.setdefault(k, []).append(model)

		"""PARTIAL Dictionary with ALL hits"""
		if partial_metamatch:
			for k, v in partial_metamatch.items():	
				for k2, v2 in partial_metamatch2.items():
					for x in v:
						partial_metamatch4.setdefault(k, []).append(x)
					for x2 in v2:
						partial_metamatch4.setdefault(k2, []).append(x2)
		elif not partial_metamatch:
			for k2, v2 in partial_metamatch2.items():
				for x2 in v2:
					partial_metamatch4.setdefault(k2, []).append(x2)						
	
		"""Creates PARTIAL dictionary without duplicate values (all hits)"""
		for k, v in partial_metamatch4.items():
			for i in range(0, len(v)):
				if v[i] not in v[i+1:]:
					partial_metamatch3.setdefault(k, []).append(v[i])

		"""PARTIAL efflux pump meta-models with snp"""
		for key, value in partial_metamatch3.items():
			partial_jsondict ={}
			x3 = {}
			x2 = {}
			bitscore_dict ={}
			percidentity_dict ={}
			partial_modellist=[]			
			for k1,v1 in components.items():
				if key == k1:
					for x in value:
						newdict = {}
						partial_modellist.append(x["model_id"])
						partial_modellist1 = set(partial_modellist)
						if all(y in partial_modellist1 for y in v1):
							if any("snp" in x for x in value):
								newdict = partial_metamatch3[key]
								count=0
								for x in newdict:
									if x["type_match"] == "Perfect" and ("model_id: "+x["model_id"]) not in x2.keys():
										bitscore_dict.setdefault((x["model_id"]), []).append(x["bit_score"])
										bitscore_dict = dict((k, list(set(v))) for k, v in bitscore_dict.items())
										percidentity_dict.setdefault((x["model_id"]), []).append(x["perc_identity"])
										percidentity_dict = dict((k, list(set(v))) for k, v in percidentity_dict.items())
										for k,v in bitscore_dict.items():
											if x["bit_score"] == max(v):
												for k,v in percidentity_dict.items():
													if x["perc_identity"] == max(v):
														x2 = {"model_id: "+x["model_id"]: x}
									elif x["type_match"] == "Strict" and ("model_id: "+x["model_id"]) not in x2.keys():
										if "snp" not in x:
											bitscore_dict.setdefault((x["model_id"]), []).append(x["bit_score"])
											bitscore_dict = dict((k, list(set(v))) for k, v in bitscore_dict.items())
											percidentity_dict.setdefault((x["model_id"]), []).append(x["perc_identity"])
											percidentity_dict = dict((k, list(set(v))) for k, v in percidentity_dict.items())
											for k,v in bitscore_dict.items():
												if x["bit_score"] == max(v):
													for k,v in percidentity_dict.items():
														if x["perc_identity"] == max(v):
															x2 = {"model_id: "+x["model_id"]: x}
										if "snp" in x:				
											if x not in x2.values():
												count+=1
												x2 = {"model_id: "+x["model_id"]+"_"+str(count): x}													
									elif x["type_match"] == "Loose" and ("model_id: "+x["model_id"]) not in x2.keys():	
										bitscore_dict.setdefault((x["model_id"]), []).append(x["bit_score"])
										bitscore_dict = dict((k, list(set(v))) for k, v in bitscore_dict.items())
										percidentity_dict.setdefault((x["model_id"]), []).append(x["perc_identity"])
										percidentity_dict = dict((k, list(set(v))) for k, v in percidentity_dict.items())
										for k,v in bitscore_dict.items():
											if x["bit_score"] == max(v):
												for k,v in percidentity_dict.items():
													if x["perc_identity"] == max(v):
														x2 = {"model_id: "+x["model_id"]: x}
									x3.update(x2)																		
									partial_jsondict["Partial"] = x3
									for entry in card_data:
										if card_data[entry]["model_id"] == key:
											if ("meta-model ID: "+key+", "+card_data[entry]["model_name"]) not in effluxpumpmetamodels.keys():
												effluxpumpmetamodels["meta-model ID: "+key+", "+card_data[entry]["model_name"]] = partial_jsondict

		"""PARTIAL efflux pump meta-models with no snp"""
		for key, value in partial_metamatch3.items():
			partial_jsondict ={}
			x3 = {}
			x2 = {}
			bitscore_dict ={}
			percidentity_dict ={}	
			perflist=[]			
			for k1,v1 in components.items():
				if key == k1:				
					for x in value:
						newdict = {}
						if x["type_match"] == "Perfect":
							perflist.append(x["model_id"])
						if (len(perflist)) >= 1:
							partial_modellist.append(x["model_id"])
							partial_modellist1 = set(partial_modellist)					
						if all(y in partial_modellist1 for y in v1):
							if all("snp" not in x for x in value):
								newdict = partial_metamatch3[key]
								for x in newdict:
									if x["type_match"] == "Perfect":									
										if ("model_id: "+x["model_id"]) not in x3.keys():
											bitscore_dict.setdefault((x["model_id"]), []).append(x["bit_score"])
											bitscore_dict = dict((k, list(set(v))) for k, v in bitscore_dict.items())
											percidentity_dict.setdefault((x["model_id"]), []).append(x["perc_identity"])
											percidentity_dict = dict((k, list(set(v))) for k, v in percidentity_dict.items())
											for k,v in bitscore_dict.items():
												if x["bit_score"] == max(v):
													for k,v in percidentity_dict.items():
														if x["perc_identity"] == max(v):
															x2 = {"model_id: "+x["model_id"]: x}	
									elif x["type_match"] == "Strict":
										if ("model_id: "+x["model_id"]) not in x3.keys():
											bitscore_dict.setdefault((x["model_id"]), []).append(x["bit_score"])
											bitscore_dict = dict((k, list(set(v))) for k, v in bitscore_dict.items())
											percidentity_dict.setdefault((x["model_id"]), []).append(x["perc_identity"])
											percidentity_dict = dict((k, list(set(v))) for k, v in percidentity_dict.items())
											for k,v in bitscore_dict.items():
												if x["bit_score"] == max(v):
													for k,v in percidentity_dict.items():
														if x["perc_identity"] == max(v):
															x2 = {"model_id: "+x["model_id"]: x}
									elif x["type_match"] == "Loose":
										if ("model_id: "+x["model_id"]) not in x3.keys():
											bitscore_dict.setdefault((x["model_id"]), []).append(x["bit_score"])
											bitscore_dict = dict((k, list(set(v))) for k, v in bitscore_dict.items())
											percidentity_dict.setdefault((x["model_id"]), []).append(x["perc_identity"])
											percidentity_dict = dict((k, list(set(v))) for k, v in percidentity_dict.items())
											for k,v in bitscore_dict.items():
												if x["bit_score"] == max(v):
													for k,v in percidentity_dict.items():
														if x["perc_identity"] == max(v):
															x2 = {"model_id: "+x["model_id"]: x}
									x3.update(x2)								
									partial_jsondict["Partial"] = x3
									for entry in card_data:
										if card_data[entry]["model_id"] == key:
											if ("meta-model ID: "+key+", "+card_data[entry]["model_name"]) not in effluxpumpmetamodels.keys():
												effluxpumpmetamodels["meta-model ID: "+key+", "+card_data[entry]["model_name"]] = partial_jsondict							

		"""PUTATIVE efflux pump meta-models"""
		
		for key, value in partial_metamatch3.items():
			putative_jsondict ={}
			x3 = {}
			x2 = {}
			bitscore_dict ={}
			percidentity_dict ={}
			for k1,v1 in components.items():
				if key == k1:
					for x in value:	
						if all("snp" not in x for x in value):
							if x["type_match"] == "Perfect":									
								if ("model_id: "+x["model_id"]) not in x3.keys():
									bitscore_dict.setdefault((x["model_id"]), []).append(x["bit_score"])
									bitscore_dict = dict((k, list(set(v))) for k, v in bitscore_dict.items())
									percidentity_dict.setdefault((x["model_id"]), []).append(x["perc_identity"])
									percidentity_dict = dict((k, list(set(v))) for k, v in percidentity_dict.items())
									for k,v in bitscore_dict.items():
										if x["bit_score"] == max(v):
											for k,v in percidentity_dict.items():
												if x["perc_identity"] == max(v):
													x2 = {"model_id: "+x["model_id"]: x}	
							elif x["type_match"] == "Strict":
								if ("model_id: "+x["model_id"]) not in x3.keys():
									bitscore_dict.setdefault((x["model_id"]), []).append(x["bit_score"])
									bitscore_dict = dict((k, list(set(v))) for k, v in bitscore_dict.items())
									percidentity_dict.setdefault((x["model_id"]), []).append(x["perc_identity"])
									percidentity_dict = dict((k, list(set(v))) for k, v in percidentity_dict.items())
									for k,v in bitscore_dict.items():
										if x["bit_score"] == max(v):
											for k,v in percidentity_dict.items():
												if x["perc_identity"] == max(v):
													x2 = {"model_id: "+x["model_id"]: x}
							elif x["type_match"] == "Loose":
								if ("model_id: "+x["model_id"]) not in x3.keys():
									bitscore_dict.setdefault((x["model_id"]), []).append(x["bit_score"])
									bitscore_dict = dict((k, list(set(v))) for k, v in bitscore_dict.items())
									percidentity_dict.setdefault((x["model_id"]), []).append(x["perc_identity"])
									percidentity_dict = dict((k, list(set(v))) for k, v in percidentity_dict.items())
									for k,v in bitscore_dict.items():
										if x["bit_score"] == max(v):
											for k,v in percidentity_dict.items():
												if x["perc_identity"] == max(v):
													x2 = {"model_id: "+x["model_id"]: x}
							x3.update(x2)								
							putative_jsondict["Putative"] = x3
							for entry in card_data:
								if card_data[entry]["model_id"] == key:
									if ("meta-model ID: "+key+", "+card_data[entry]["model_name"]) not in effluxpumpmetamodels.keys():
										effluxpumpmetamodels["meta-model ID: "+key+", "+card_data[entry]["model_name"]] = putative_jsondict




		ejson = json.dumps(effluxpumpmetamodels)		
		newfilename = os.path.splitext(newrgi_filename)[0]
		if effluxpumpmetamodels:
			with open(filepath+"/"+"EffluxPumpMetaModelResults"+'/'+newfilename+'.'+'effluxpumpmetamodels.json', 'w') as f:
				f.write(ejson)

		for k,v in effluxpumpmetamodels.items():
			modellist =[]
			for k1,v1 in v.items():
				for k2,v2 in v1.items():
					modellist.append(v2["model_id"])
#			logger.info(k+" | "+k1+" | model ID:"+str((", ".join(modellist))))	

filepath = os.path.dirname(os.path.realpath(__file__))
if not os.path.exists(filepath+"/"+"EffluxPumpMetaModelResults"):
	os.makedirs("EffluxPumpMetaModelResults")

for file in os.listdir(filepath):
	if file.endswith(".json"):
		obj=EffluxPump(os.path.join(file))
		obj.run()
