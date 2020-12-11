from Bio import AlignIO
from Bio import Phylo
from Bio.Align import substitution_matrices
from Bio.Align.Applications import ClustalwCommandline
from io import StringIO
from numpy import NaN
import pandas as pd
import re


# Node class used in UPGMA algorithm
class UPGMANode:
	def __init__(self, name):
		self.name = name
		self.age = 0
		self.children = []
		self.parent = None

	def addChild(self, child):
		self.children.append(child)

	def setParent(self, parent):
		self.parent = parent

	def setAge(self, age):
		self.age = age

	def countDescendants(self):
		if len(self.children) == 0:
			return 1
		
		count = 0
		for child in self.children:
			count += child.countDescendants()
		return count

	def toStringH(self, previous):
		s = ""
		if self.parent is not None and self.parent is not previous:
			s += self.parent.toStringH(self) + ":" + str(abs(self.parent.age - self.age)) + ", "
		for child in self.children:
			if child is not previous:
				s += child.toStringH(self) + ":" + str(abs(child.age - self.age)) + ", "

		s = "(" + s[:-2] + ")" #there's probably a better way to do this but I don't care
		s = s.replace("()","")
		if isinstance(self.name, str): #this is so we can ignore intermidiary nodes when drawing phylogenetic trees
			s += self.name
		
		return s

	def toString(self):
		return self.toStringH(None)


#UPGMA algorithm for generating a phylogenetic tree from the provided distance matrix
def UPGMA(distances):
	isDataFrame = isinstance(distances, pd.DataFrame)
	labels = distances.columns if isDataFrame else range(len(distances))

	if not isDataFrame:
		distances = pd.DataFrame(
			data = distances,
			index = labels,
			columns = labels)
	
	for cluster in distances.columns: #We want to ignore the diagonal, so to make sure no math is done with it we set the diagonal to NaN
		distances[cluster][cluster] = NaN

	nodes = {}

	for name in labels:
		nodes[name] = UPGMANode(name)

	Cnew = len(distances)

	#initialization done
	
	while distances.size > 1:
		Ci, Cj = distances.stack().idxmin() #find the two closest clusters Ci abd Cj (these are indices, not the clusters themselves) 

		distances[Cnew] = ((distances[Ci] * nodes[Ci].countDescendants()) + (distances[Cj] * nodes[Cj].countDescendants())) / (nodes[Ci].countDescendants() + nodes[Cj].countDescendants()) #merge Ci and Cj into a new cluster Cnew with |Ci|+|Cj| elements
		distances.loc[Cnew] = distances[Cnew]

		nodes[Cnew] = UPGMANode(Cnew) #add a new node labeled by cluster Cnew to T (nodes)

		nodes[Cnew].addChild(nodes[Ci]) #connect node Cnew to Ci and Cj by directed edges
		nodes[Cnew].addChild(nodes[Cj])
		nodes[Ci].setParent(nodes[Cnew])
		nodes[Cj].setParent(nodes[Cnew])

		nodes[Cnew].setAge(distances[Ci][Cj]) #Age(Cnew) <- D[C1,C2] (supposed to divide by 2 but there's no mathematical reason two. The relative proportions are the same)

		distances = distances.drop(Ci, axis=0) #remove the rows and columns of D corresponding to Ci and Cj
		distances = distances.drop(Cj, axis=0)
		distances = distances.drop(Ci, axis=1)
		distances = distances.drop(Cj, axis=1)

		Cnew += 1

	s = nodes[distances.columns[0]].toString()
	tree = Phylo.read(StringIO(s), "newick")
	Phylo.draw(tree)

	#s = "" #CODE FOR ROSALIND OUTPUT
	#for node in list(nodes.values()): #  for each edge (v, w) in T (nodes) length of (v, w) <- Age(v) - Age(w) 
	#	if node.parent is not None:
	#		s += str(node.name) + "->" + str(node.parent.name) + ":" + str(node.parent.age - node.age) + "\n"
	#	for child in node.children:
	#		s += str(node.name) + "->" + str(child.name) + ":" + str(abs(child.age - node.age)) + "\n"
	#return s[:-1]


# Build a distance matrix using BLOSUM from the provided multiple sequence alignment
# NOTE: this uses the NEGATIVE of the BLOSUM scores so it can be used with existing UPGMA algorithm
def getBLOSUMDistanceMatrix(alignment):
	blosumMatrix = substitution_matrices.load("BLOSUM62")
	
	df = pd.DataFrame(columns=list(r.id for r in alignment),index=list(r.id for r in alignment))

	for record1 in alignment:
		for record2 in alignment:
			score = 0
			for i in range(len(record1.seq)):
				aa1 = record1[i] if record1[i] != '-' else '*'
				aa2 = record2[i] if record2[i] != '-' else '*'
				score -= blosumMatrix[aa1][aa2]
			df[record1.id][record2.id] = score

	return df.apply(pd.to_numeric)
	

# create a multiple sequence alignment from provided file
def msa(filePath):
	ClustalwCommandline("clustalw2", infile=filePath, clustering="UPGMA")()
	return AlignIO.read(filePath.split(".")[0] + ".aln", "clustal")


#rename fasta strings in file for species names (this is a QOL function used for testing)
def renameFasta(inputFilePath, outputFilePath):
	newFile = open(outputFilePath, "w")
	s = ""
	for line in open(inputFilePath, "r").readlines():
		start = line.find('[')
		end = line.find(']')
		s += ">"+line[start+1:end].replace(" ", "_")+"\n" if start > -1 else line
	newFile.write(s)

# test code
#buildTree(UPGMA(fileToMatrix("rosalind_ba7d.txt")))
#buildTree(UPGMA(fileToMatrix("UPGMA.txt")))

renameFasta("h1raw.fasta", "h1.fasta")
UPGMA(getBLOSUMDistanceMatrix(msa("h1.fasta")))
tree = Phylo.read("h1.dnd", "newick")
Phylo.draw(tree)

#renameFasta("h4raw.fasta", "h4.fasta")
#UPGMA(getBLOSUMDistanceMatrix(msa("h4.fasta")))

#renameFasta("megaTestraw.fasta", "megaTest.fasta")
#UPGMA(getBLOSUMDistanceMatrix(msa("megaTest.fasta")))
#tree = Phylo.read("megaTest.dnd", "newick")
#Phylo.draw(tree)





#================================================== Unused code kept for posterity ==================================================


#class PhyloNode:
#	def __init__(self, name):
#		self.name = name
#		self.children = {}
#
#
#	def addChild(self, node, weight):
#		self.children[node] = weight
#
#
#	def toStringH(self, parent):
#		if len(self.children) == 1:
#			if parent in self.children:
#				return self.name
#			else:
#				return list(self.children.keys())[0].toStringH(self) + ":" + str(list(self.children.values())[0])
#		
#		s = "("
#		for child in self.children:
#			if child is not parent:
#				s = s + child.toStringH(self) + ":" + str(self.children[child]) + ", "
#		
#		return s[:-2] + ")" + self.name
#
#
#	def toString(self):
#		if len(self.children) == 1:
#			return "(" + self.toStringH(None) + ")" + self.name
#		
#		return self.toStringH(self)
#		
#
#def buildTree(graphString):
#	treeDict = {} #num: node
#
#	edges = re.split("\n", graphString)
#	
#	for edge in edges:
#		vals = re.split("->|:", edge)
#				
#		node1ID = vals[0]
#		node2ID = vals[1]
#		weight = vals[2]
#
#		if node1ID not in treeDict:
#			treeDict[node1ID] = PhyloNode(node1ID)
#
#		if node2ID not in treeDict:
#			treeDict[node2ID] = PhyloNode(node2ID)
#
#		treeDict[node1ID].addChild(treeDict[node2ID], weight)
#
#	s = list(treeDict.values())[0].toString()
#	tree = Phylo.read(StringIO(s), "newick")
#	Phylo.draw(tree)
#
#	
#def parseRos(filePath): #this code isn't used anymore (it was a workaround to get the rosalind UPGMA output to work)
#	treeDict = {} #num: node
#
#	for line in open(filePath, "r").readlines():
#		vals = re.findall("[A-Z]|[0-9]+", line)
#		
#		node1ID = vals[0]
#		node2ID = vals[1]
#		weight = vals[2]
#
#		if node1ID not in treeDict:
#			treeDict[node1ID] = PhyloNode(node1ID)
#
#		if node2ID not in treeDict:
#			treeDict[node2ID] = PhyloNode(node2ID)
#
#		treeDict[node1ID].addChild(treeDict[node2ID], weight)
#
#	s = treeDict['0'].toString()
#	tree = Phylo.read(StringIO(s), "newick")
#	Phylo.draw(tree)
#
#
#def fileToMatrix(filePath):
#	matrix = []
#	for line in open(filePath, "r").readlines():
#		row = []
#		for val in re.findall("[0-9]+", line):
#			row.append(int(val))
#		matrix.append(row)
#
#	return matrix
#
##parseRos("test.txt")
#
