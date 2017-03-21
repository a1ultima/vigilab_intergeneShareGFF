
import pdb

"""

DEPRECATED: 

	See: "@vigilab_intergeneGFF" tool (i.e. "alternative 2") for latest one

DESCRIPTION:

	... GFF file modifier for Kathrin. It takes ... 

"""


### alternative 1 {{

#
# 1. Read in the file
#


#
# 2. Ensure each row is ordered
#

# ... i.e. such that genes to the far left are on the top row:
# ===1====>  	<==2==   <==3==  ==4==>  
#
# row 1: ===1===> 
# row 2: <==2== 
# row 3: <==3==
# row 4: ==4==>
#

#
# 3. Start from 0 (intergene_1_start=0), and keep counting +1 until we hit the first gene's START (intergene_1_end = gene1_START ... +1 ...)
#

#
# 4. Keep counting until we reach the first gene's END: 
#

# ... intergene_1_end = gene1_END ... +1 ...

# but also store the gene1_END, to later partition the intergenic region equally between adjacent "tug of war" genes

#
# 5. Keep on counting until we finally reach the second gene's start

# ... intergene_1_end = gene2_START
#

# }} 1 alternative 2 {{


""" PROTO-README (@TODO: after testing, incorporate into ./README.md)

# TITLE: 
	@@vigilab_intergeneShareGFF tool (for Kathrin)

# DESCRIPTION:
	for Kathrin to modify a GFF input file (e.g. ./toy.gff) to output a file whose "transcript"-only rows have their "start" and "end" fields extended to include neighbouring intergenic regions (see: README.md for visual intuition)
	- So far we assume that the "types" we are interested in are "transcript" (filtering those rows accordingly) // @TODO: change later to something Kathrin wants: e.g. "cds"?
	- We also assume that there is only a single chromosome (for Kathrin there will be >1) // @TODO: later change to aknowledge >1 chromosomes

** Nomenclature: **
	# @@gff: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md,
	# @TODO: is the .gff the correct version? (gff3 or gff2?)
	# @DONE: make sure the #keys == #transcript rows
	# @@GFF_Obj: "transcript" Filtered GFF object

** Pipeline: **

	1. Read in the gff file  (@@read-gff)
		- @DONE: see what I did for the tug-o-war and pulley-seq // @A: not relevant
		- @DONE: Parse .gff file, filter according to "transcript" rows only (see: field: "type"), @TODO: later filter according to @kathrin's needs, e.g. "cds"?

	2. Filter to keep only transcripts (@@filter-gff)
		- @TODO: in Kathrin's case, it may need to be "cds" and not "transcript", change later
		- @TODO: also in Kathrin's case, we need to be aware of the chromosome of each "cds"

	3. Create gene objects:  (@@gene_objects, @@gff_obj, @@gff_i_obj)
		- e.g. a = to share sequences between themselves (e.g. a.me) and their leftwards neighbour (e.g. a.left)	
		- We want to add two new fields to the gene_to_field dict, so that each gene has a "reg_start", and a "reg_end"
			
	4. FOR loop to iterature a.share_neighbouring_seqs for all genes: (@@share-neighbours)
		- @TODO: shall we use a copy.deepcopy() in each iteration of the gene object contruction? I hope not, check memory usage

	5. Write each gene's updated "start" and "end" fields to a new gff file (@@write-gff)
		- for each a.me, write the modified "start" and "end" fields to a file (apropos a.me_new["start_intergene"], a.me_new["end_intergene"]), 
		- ..rather than doing this in the gene objects (a.me) we can create a separate function to write each instance of a.me's data to the file in one go (i.e. we do not want to make a separate method for writing to file for each a. class)
"""

#
# 1. Read in the gff file  (@@read-gff)
#

#
# 2. Filter to keep only transcripts (@@filter-gff)
#

# @TODO: wrap into a function(), [e.g. def readGffFile()] for now, here's the doctring: 

""" readGffFile (@@readGffFile, @read-gff, @filter-gff)

Description:
	- reads a .gff input file, then filters datarows so only those matching the "feature_type" arg are kept.
	- e.g. ./toy.gff
Args:
	- feature_type (str) -- filters incoming .gff file to keep only rows whose field "type" matches feature_type, e.g."transcript" (default)


"""


# gene_to_field is where all the datarows/cols in the input gff are stored, 

gene_to_field = {}  # keys: genes represented as 1..n, values: the 8 fields (cols) of a gff row

gene_i = 0

# Reading file buffer...

with open("./toy.gff", "r") as fi:

	print("Reading GFF file into: gene_to_field (dict), index as such: gene_to_field[gene_i], where gene_i is between 1-to-n...")
	
	while True:

		line = fi.readline().rstrip()

		if line == "":
			break
		
		line_split = line.split("\t")

		if line_split[2] != "transcript":
			continue

		#c1_reference_seq = line_split[0] # e.g. 'scaffold_150' 
		#c2_source = line_split[1] # e.g. 'GWSUNI'
		#c3_type = line_split[2] # e.g. 'transcript'
		#c4_start = line_split[3] # e.g. '1372'
		#c5_end = line_split[4] # e.g. '2031'
		#c6_score = line_split[5] # e.g. '45.89'
		#c7_strand = line_split[6] # e.g. '+'
		#c8_phase = line_split[7] # e.g. '.' @Note: codon frame (0,1,2)
		#c9_attributes = line_split[8] # e.g. <see @gff3.md>

		gene_i += 1  # indexing starts from 1, i.e. [1] = first gene 

		##@TEST: sometimes 4.00 instead of 4.0 (trivial)
		#if not (str(line_split[5])==str(float(line_split[5]))):
		#	print("oops")
		#	print("\t"+str(line_split[5])+"___"+str(float(line_split[5])))

		gene_to_field[gene_i] = { \
			"c1_reference_seq":line_split[0],# e.g. 'scaffold_150' \
			"c2_source":line_split[1],# e.g. 'GWSUNI' \
			"c3_type":line_split[2],# e.g. 'transcript' \
			"c4_start":int(line_split[3]),# e.g. '1372' \
			"c5_end":int(line_split[4]),# e.g. '2031' \
			"c6_score":float(line_split[5]),# e.g. '45.89' \
			"c7_strand":line_split[6],# e.g. '+' \
			"c8_phase":line_split[7],# e.g. '.' @Note: codon frame (0,1,2) \
			"c9_attributes":line_split[8]# e.g. <see @gff3.md> \
		}	
		



## @TEST:@DONE: print the contents of gene_to_field back out in gff format

#
# 2. <add the pipeline step here @TODO:>
#

###########
# Classes #
###########

# Note: I could have done this in lower-level langauge, but just in case Kathrin wanted to do other things with it in future I wanted to wrap separate tasks in modular object form.

#
# Exceptions & Errors (error handling)
#

class IntergenicSeqDiffImpossible(Exception):	
	""" Base class for exceptions involving negative (-ve) values encountered for..
		..intergenic_seq_diff variable, we never expect -ve values.

	"""
	# @LATEST-2017-03-06:@2033 - @TODO: throw this error class when -ve intergenic_seq_diff occurs in the gene_and_neighbours() class, @TODO: make sure the gene_and_neighbours() class is later consistently capitalised, PEP8	
	pass

class DuplicateFeatureError(IntergenicSeqDiffImpossible):
	""" Error raised when we think the InterGenicSeqDiffImpossible Exception flags a duplicated feature.
		
		Attributes:
			expression -- input expression in which the error occurred
			message -- explanation of the error
	"""
	
	# def __init__(self, expression, message):  # @VANILLA:@TODO:later maybe also allow for expression and message, but need to think more about how to work it 
		#self.expression = expression
		#self.message = message
	def __init__(self, message):  # @VANILLA:@TODO:later maybe also allow for expression, but need to think more about how to work it 
		# @TODO: replace ^ and ^^ (above two lines) to have the following error message instead
		self.message = message # e.g. something like: "\t"*5+"Me (gene_i=%r) and Left neighbour (gene_i=%r) have a -ve tandem difference in integenic location in BPs (intergenic_seq_diff=%r)... we may be duplicates! Writing to file... \n" % (self.my_id, self.left_id, intergenic_seq_diff)
		# @DONE: add some error logs to "./log.txt" 
		with open("./log.txt", "a") as fo_log:
			fo_log.write(self.message)

#
# Error handlers 
#
def assertNoDuplicateFeatures(intergenic_seq_diff, gffFeatureObj):
	""" Error Handling Function for checking whether or not an encounter with -ve intergenic_seq_diff values is caused specifically by duplicate features..
		..e.g. "transcript", where the "start" and "end" fields are equal.

	ARG:
		- gffFeatureObj (class) -- e.g. a = gene_and_neighbours() class instance

	"""
	
	gene_i_me = gffFeatureObj.my_id # e.g. 8 (int), in: gene_and_neighbours(8)
	gene_i_left = gffFeatureObj.left_id # e.g. 7 (int), ^ since: gene_i_left = gene_i_me - 1

	if intergenic_seq_diff<0:  # -ve values of intergenic_seq_diff are impossible	@LATEST:@2221
		if ((gffFeatureObj.me["c4_start"]==gffFeatureObj.left["c4_start"]) and (gffFeatureObj.me["c5_end"]==gffFeatureObj.left["c5_end"])): # the "start" and "end" fields of me and left are equal, we are the same (I hope) @TODO:more throrough testing by comparing equality of all fields)				
			raise DuplicateFeatureError("\t"*5+"Me (gene_i=%r) and Left neighbour (gene_i=%r) have a -ve tandem difference in integenic location in BPs (intergenic_seq_diff=%r)... we may be duplicates! Writing to file... \n" % (gene_i_me, gene_i_left, intergenic_seq_diff))
 # @LATEST:2017-03-06:@2212:@TODO:consistent naming
		else:
			raise IntergenicSeqDiffImpossible  
	else:
		print("\t"*5+"No problems so far...")  # @TODO: so we've dodged the bug
#
# Objects
#

class gene_and_neighbours(object):
	
	""" @TODO: rename to: "@@gffFeatureObj"
	
	Description:
		- Object for handling incoming gff input file data, ..
		..representing a row of data whose type is specified as "transcript" @TODO:change so that type can be specified through some input argument to the class instantiation, e.g. "self.feature_type". 
	
	Args:
		- self.gene_i (int) -- number that indexes feature_type-filtered data row of: @gene_to_field, as in gene_to_field[gene_i], e.g. gene_i = 2, where 2 is te 2nd datarow of gene_to_field


 	""" 


	# @TODO:run a test to ensure the class of the gene_i_to_fields is an indexed @GFF_OBJ, e.g. gene_to_fields[i], where i enumerate 0:len(gene_to_fields)
	
	def __init__(self, gene_i):  # see: @GFF_OBJ (@TODO: make the "see:" description point to file if this class file is in a separate file)
		
		""" Parse the field data attributed to a single gene, and also its left/right neighbours (@TODO: do we keep this in the same file as the gene_to_field parser?) """

		print("\tGene_i: "+str(gene_i))
		print("\t\tCreating myself (a.me) and my leftward neighbour (a.left)...")
		
		#
		#       [     -1     ] [     0      ] [      +1     ] 
		# e.g.  <--self.left-- --self.self--> --self.right-->
		#
		self.me = gene_to_field[gene_i]    # Parse the field data attributed to the neighbouring gene (whose "start" < self."end") 

		self.my_id = gene_i
		self.left_id = gene_i-1
		
		# @TODO: catch the edge case, where gene 1 has no left neighbour, and gene n has no right neighbour
		# @TODO: deal with the cases where there are multiple chromosomes! 

		if (gene_i-1) == 0:
			print("Gene to the left of gene_i="+str(gene_i)+" does not exist, since gene_i is the furthest left in the chromosome, setting to None...")
			self.left = None
		else:
			self.left = gene_to_field[gene_i-1]
		
		## @DONE: no need to parse the right-neighbour, remove
		#if (gene_i+1) == (len(gene_to_field.keys())+1):
		#	print("Gene to the right of gene_i="+str(gene_i)+" does not exist, since gene_i is the furthest right in the chromosome, setting to None")
		#	self.right = None
		#else: 
		#	self.right = gene_to_field[gene_i+1]	

		print("\t\t\tSelf data successfully loaded from gene_to_field (dict) [into self.me]...")
	
	def share_neighbouring_seqs(self):

		""" Determine which orientations the neighbouring seqs are in: head-to-head? head-to-tail? tail-to-head? tail-to-tail? """

		print("\t\tAttempting to share sequences between me (a.me) and my leftward neighbour...")

		if self.left==None:
			# @TODO: make sure to deal with left-edge case
			pdb.set_trace()
	
		# @Note:@@factored-out-@intergenic_seq_diff: two reasons: (i) catch intergenic_seq_diff == -ve value errors, and (ii) reduce code redundancy, 
		intergenic_seq_diff = self.me["c4_start"] - self.left["c5_end"]  # @TODO: code repetition, can factor out
		print("\t\t\tTotal No. intergenic BPs between me's START field and left's END field: "+str(intergenic_seq_diff))

		# @DEBUG:@TODO: If the a.left neighbour is a duplicate of a.me we get
		# .. a -ve value. 
		# 	- @Q: How can we deal with these duplicates best?
		# 		- @UNDONE:@A: Used an assertion statement that breaks the program for now
		# 		- @DONE:@A: @TODO: need to catch it and use as warning instead? 
		# 	 // @LATEST:2017-03-06-@1929

		# @DuplicateFeatureError handling, to test if me and left neighbour are the same feature (duplicated)
		print("\t"*4+"Checking for bugs...")
		try:
			assertNoDuplicateFeatures(intergenic_seq_diff, self) # error logs written to: ./log.txt
		except DuplicateFeatureError as err:
			print(err.message) 
			print("\t"*5+"...Skipping gene_i: %r" % gene_i)
			pass  # @TODO: better to skip it // instead of skipping, it currently: carries on processing as if no err was encountered

		# @DONE: @LATEST-2017-03-06-1900: we factored out the intergenic_seq_diff

		#############################################
		# CASE 1 & 2: Tail-to-Tail and Head-to-Head #
		#############################################

		# Share seqs with left, Head-to-head (1/2 each), tail-to-tail (2/3 each)
		# @TODO: can reduce code by doing an IF strand directions are not the same then share 50:50 ...
		if ((self.left["c7_strand"]=="+") and (self.me["c7_strand"]=="-")) or ((self.left["c7_strand"]=="-") and (self.right["c7_strand"]=="+")):
			print("\t\tHead-to-Head (or Tail-to-Tail) case encountered: 5'===left===>3'......3'<===me===5'")	
			# @TODO: I assume the gff convention is to name the left-most chromosome pos as 0

			## @NOTE: factored out into upper indent: see: @@factored-out-@intergenic_seq_diff
			#intergenic_seq_diff = self.me["c4_start"] - self.left["c5_end"]  # @TODO: code repetition, can factor out
			#print("\t\t\tNo. intergenic BPs between me and left: "+str(intergenic_seq_diff))
			
			fraction_me = (1/2)
			print("\t\t\tFraction of BPs taken by me: "+str(fraction_me)) 
			
			me_share   = intergenic_seq_diff * fraction_me
			print("\t\t\tMy share of BPs: "+str(me_share)) 
			
			left_share = intergenic_seq_diff * (1-fraction_me)
			print("\t\t\tLeft Neighbour share of BPs: "+str(left_share)) 

			# @LATEST:@TODO:add the shared intergenic positions to the self.me, or directly to file?
			
			# @TEST:@DONE: it ^, done: indeed correct orientation found and also correct sharing fractions

		########################
		# CASE 3: Tail-to-head #
		########################

		if ((self.left["c7_strand"]=="+") and (self.me["c7_strand"]=="+")):
			print("\t\tTail-to-Head case encountered: 5'===left===>3'...intergenic...5'===me===>3'")

			pdb.set_trace() # GWSUNIT00000016001 // @ANDY-and-Luke: python tutorial: how to debug // @Luke PDBTutorial:@DATE:@2017-03-20

			## @NOTE: factored out into upper indent: see: @@factored-out-@intergenic_seq_diff
			#intergenic_seq_diff = self.me["c4_start"] - self.left["c5_end"]  # @TODO: code repetition, can factor out
			#print("\t\t\tNo. intergenic BPs between me and left: "+str(intergenic_seq_diff))

			fraction_me = (2/3) # @DONE: check with Kathrin: 2/3 to 5' or 2/3 to the 3'?	
			print("\t\t\tFraction of BPs taken by me: "+str(fraction_me)) 

			me_share   = intergenic_seq_diff * fraction_me
			print("\t\t\tMy share of BPs: "+str(me_share)) 

			left_share = intergenic_seq_diff * (1-fraction_me)
			print("\t\t\tLeft Neighbour share of BPs: "+str(left_share)) 

			# @TEST: it^i  


		########################
		# CASE 4: Tail-to-head #
		########################

		if ((self.left["c7_strand"]=="-") and (self.me["c7_strand"]=="-")):
			print("\t\tHead-to-Tail case encountered: 3'<===left===5'...intergenic...3'<===me===5'")

			## @NOTE: factored out into upper indent: see: @@factored-out-@intergenic_seq_diff
			#intergenic_seq_diff = self.me["c4_start"] - self.left["c5_end"]  # @TODO: code repetition, can factor out
			#print("\t\t\tNo. intergenic BPs between me and left: "+str(intergenic_seq_diff))
			
			fraction_me = (1/3)
			print("\t\t\tFraction of BPs taken by me: "+str(fraction_me)) 

			me_share = intergeneic_seq_diff * fraction_me 
			print("\t\t\tMy share of BPs: "+str(me_share)) 

			left_share = intergenic_seq_diff * (1-fraction_me)
			print("\t\t\tLeft Neighbour share of BPs: "+str(left_share))

			pdb.set_trace() # @2227 // @TODO:REMOVE 

			# @TEST: it^ // @LATEST:2017-03-20-@2225

#
# 3. Create gene objects: (see: @gene_objects, @gff_obj, @gff_i_obj)
#


## @DONE:@TEST: example case: Tail-to-Tail: 5'(head)===>3'(tail)...|...3'(tail)<===5'(head)..
# ...All numbers consistent to manually calculated example (./toy.gff)
#a = gene_and_neighbours( 2 )
#a.share_neighbouring_seqs()

## @TODO:@TEST: example case: Tail-to-Head: 5'(head)===>3'(tail)..|....5'(head)===>3'(tail)..

## @NOTES:
# 	-@TODO: catch the error, throw it to use as msg, then skip the error-causing case // @DONE: catch errors, and just "pass" as if error was not there // There are multiple duplicate transcripts, these cause intergeneic_seq_diff to be -ve, we need to catch -ve ones and throw an error

# ...@LATEST-2017-03-06.. @LATEST-2017-03-20

b = gene_and_neighbours( 8 )  # @Q: what does this do // @A: see: gene_and_share

b.share_neighbouring_seqs()


# }} 2 alternative 

