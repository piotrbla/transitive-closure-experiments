values = []
i = 2
with open('results.txt') as f:
	for line in f:
		if line.startswith('1600;'):
				s = line[5:].rstrip().replace(',', '.')
				values.append((i, float(s)))
				i += 2
for i, v in values:
	if (v<4.32):
		print (i, v)