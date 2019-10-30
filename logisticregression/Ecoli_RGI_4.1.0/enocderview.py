import pickle
with open('rgi_encoded.pkl', 'rb') as f:
	data = pickle.load(f)
	print(data)