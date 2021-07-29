from sklearn.tree import DecisionTreeClassifier
from sklearn.tree import export_text
import joblib
import sys
from sklearn.tree import _tree

sys.setrecursionlimit(4000)

loaded_model = joblib.load(sys.argv[1])
headers = joblib.load(sys.argv[2])
output_file = sys.argv[3]

finalHeaders = []

with open(output_file, 'r') as f:
	for line in f:
		if line[0] == "[" and len(line) > 1000:
			print(line)
			finalHeaders = line[1:-1].split(",")

			for i in range(len(finalHeaders)):
				finalHeaders[i] = finalHeaders[i].strip()
				finalHeaders[i] = finalHeaders[i].replace("'", "")

print(finalHeaders)

categories = ['A', 'C', 'G', 'T', 'N', '-']

classes = loaded_model.classes_


def getClass(weightArray):
	weightArray = weightArray[0]

	maxValue = 0
	maxIndex = 0

	for i in range(len(weightArray)):
		if weightArray[i] > maxValue:
			maxValue = weightArray[i]
			maxIndex = i

	return classes[maxIndex]


def reformatRule(sign, name, threshold):
	position = name.split("_")[0]
	letter = name.split("_")[1]

	symbol = "=="
	if sign == "<=":
		symbol = "!="

	return position + symbol + "'" + letter + "'"


outcomesToRules = dict()


def tree_to_code(tree, feature_names):
    tree_ = tree.tree_
    feature_name = [
        feature_names[i] if i != _tree.TREE_UNDEFINED else "undefined!"
        for i in tree_.feature
    ]

    prevRules = []

    def recurse(node, depth, prevRules):
        indent = "  " * depth
        if tree_.feature[node] != _tree.TREE_UNDEFINED:
            name = feature_name[node]
            threshold = tree_.threshold[node]

            prevRulesRight = prevRules.copy()
            prevRulesLeft = prevRules.copy()

            prevRulesLeft.append(reformatRule("<=", name, threshold))
            prevRulesRight.append(reformatRule(">", name, threshold))

            recurse(tree_.children_left[node], depth + 1, prevRulesLeft)
            recurse(tree_.children_right[node], depth + 1, prevRulesRight)
        else:
            rulesString = prevRules[0]
            for rule in prevRules[1:]:
            	rule = rule.strip()

            	rulesString = rulesString + "," + rule

            outcome = getClass(tree_.value[node])

            if outcome not in outcomesToRules:
            	outcomesToRules[outcome] = []

            print(outcome + "\t" + rulesString)

    recurse(0, 1, prevRules)

tree_to_code(loaded_model, finalHeaders)


def formatCommonRules(ruleSets):
	splitRuleSets = []
	totalRules = []
	commonRules = []

	for r in ruleSets:
		splitRuleSets.append(r.split(","))

		for i in r.split(","):
			totalRules.append(i)

	for i in totalRules:

		addRule = True
		for r in splitRuleSets:
			if i not in r:
				addRule = False

		if addRule:
			commonRules.append(i)

			for r in splitRuleSets:
				r.remove(i)

	retString = ""

	if len(commonRules) > 0:
		retString = commonRules[0]

		for c in commonRules[1:]:
			retString = retString + "," + c

		retString = retString + " AND "

	for s in splitRuleSets:
		if len(s) > 0:
			retString = retString + "(" + s[0]
			for r in s[1:]:
				retString = retString + "," + r

			retString = retString + ") OR "


	return retString[:-3]
