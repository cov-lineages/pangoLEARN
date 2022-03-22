import sys

outputfile = sys.argv[1]

nameToLineage = dict()

class Lineage:
	def __init__(self, name):
		self.name = name
		self.precisions = []
		self.recalls = []
		self.f1s = []
		self.supports = []

	def printStats(self):
		avgPrec = sum(self.precisions) / len(self.precisions)
		avgRec = sum(self.recalls) / len(self.recalls)
		avgF1 = sum(self.f1s) / len(self.f1s) 
		totalSupport = sum(self.supports)

		print(self.name + "," + str(avgPrec) + "," + str(avgRec) + "," + str(avgF1) + "," + str(totalSupport))

with open(outputfile, 'r') as f:
	for line in f:
		line = line.strip()
		split = line.split()

		if "macro avg" in line or "weighted avg" in line:

			name = split[0] + " " + split[1]

			if name not in nameToLineage:
				nameToLineage[name] = Lineage(name)

			nameToLineage[name].precisions.append(float(split[2]))
			nameToLineage[name].recalls.append(float(split[3]))
			nameToLineage[name].f1s.append(float(split[4]))
			nameToLineage[name].supports.append(int(split[5]))

		if len(split) == 5 and ":" not in line and "read" not in line:
			name = split[0]

			if name not in nameToLineage:
				nameToLineage[name] = Lineage(name)

			nameToLineage[name].precisions.append(float(split[1]))
			nameToLineage[name].recalls.append(float(split[2]))
			nameToLineage[name].f1s.append(float(split[3]))
			nameToLineage[name].supports.append(int(split[4]))

print("lineage,precision,recall,f1_score,support")

for key in nameToLineage:
	if "macro" not in key and "weighted" not in key:
		nameToLineage[key].printStats()

nameToLineage["macro avg"].printStats()
nameToLineage["weighted avg"].printStats()
