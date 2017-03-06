
import pdb

"""

DEPRECATED: 

	See: "vigilab_intergeneGFF" tool (i.e. "alternative 2") for new

DESCRIPTION:

	Kathrin's GFF file modifier. It takes ... 

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


""" Kathrin's vigilab_intergeneShareGFF tool for modifying a GFF input file (e.g. ./toy.gff) to output a file whose "transcript"-only rows have their "start" and "end" fields extended to include neighbouring intergenic regions (see: README.md for visual intuition)

@DESCRIPTION:
	- So far we assume that the "types" we are interested in are "transcript" (filtering those rows accordingly)
	- We also assume that there is only a single chromosome (for Kathrin there will be >1)

Nomenclature:
	# @@gff: https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md,
	# @TODO: is the .gff the correct version? (gff3 or gff2?)
	# @DONE: make sure the #keys == #transcript rows
	# @@GFF_Obj: "transcript" Filtered GFF object

Pipeline:

	# 1. Read in the gff file  (@@read-gff)
		# @DONE: see what I did for the tug-o-war and pulley-seq // @A: not relevant

	# 2. Filter to keep only transcripts (@@filter-gff)
		# @TODO: in Kathrin's case, it may need to be "cds" and not "transcript", change later
		# @TODO: also in Kathrin's case, we need to be aware of the chromosome of each "cds"

	# 3. Create gene objects:  (@@gene_objects, @@gff_obj, @@gff_i_obj)
		# e.g. a = to share sequences between themselves (e.g. a.me) and their leftwards neighbour (e.g. a.left)	

	# 4. FOR loop to iterature a.share_neighbouring_seqs for all genes: (@@share-neighbours)

		# @TODO: shall we use a copy.deepcopy() in each iteration of the gene object contruction? I hope not, check memory usage

	# 5. Write each gene's updated "start" and "end" fields to a new gff file

		# for each a.me, write the modified "start" and "end" fields to a file (apropos a.me_new["start_intergene"], a.me_new["end_intergene"]), 

		# ..rather than doing this in the gene objects (a.me) we can create a separate function to write each instance of a.me's data to the file in one go (i.e. we do not want to make a separate method for writing to file for each a. class)
"""

#
# Parse .gff file, filter according to "transcript" rows only (see: field: "type"), @TODO: later filter according to @kathrin's needs, e.g. "cds"?
#

# gene_to_field is where all the datarows/cols in the input gff are stored, 

gene_to_field = {}  # keys: genes represented as 1..n, values: the 8 fields (cols) of a gff row

gene_i = 0

#
# Reading file buffer 
#

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
		

## @TEST: print the contents of gene_to_field back out in gff format

# Note: I could have done this in lower-level langauge, but just in case kathrin wanted to do other things with it in future I wanted to wrap separate tasks in modular object form.

class gene_and_neighbours(object):

# @TODO:run a test to ensure the class of the gene_i_to_fields is an indexed @GFF_OBJ, e.g. gene_to_fields[i], where i enumerate 0:len(gene_to_fields)
	def __init__(self, gene_i):  # see: @GFF_OBJ (@TODO: make the "see:" description point to file if this class file is in a separate file)
		
		""" Parse the field data attributed to a single gene, and also its left/right neighbours (@TODO: do we keep this in the same file as the gene_to_field parser?) """

		print("\tCreating myself and my leftward neighbour...")
		
		#
		#       [     -1     ] [     0      ] [      +1     ] 
		# e.g.  <--self.left-- --self.self--> --self.right-->
		#
		self.me = gene_to_field[gene_i]    # Parse the field data attributed to the neighbouring gene (whose "start" < self."end") 

		print("\t\tSelf data loaded from gene_to_field (dict) [into self.me]...")

		# @TODO: catch the edge case, where gene 1 has no left neighbour, and gene n has no right neighbour
		# @TODO: deal with the cases where there are multiple chromosomes! 

		if (gene_i-1) == 0:
			print("Gene to the left of gene_i="+str(gene_i)+" does not exist, since gene_i is the furthest left in the chromosome, setting to None...")
			self.left = None
		else:
			self.left = gene_to_field[gene_i-1]

		if (gene_i+1) == (len(gene_to_field.keys())+1):
			print("Gene to the right of gene_i="+str(gene_i)+" does not exist, since gene_i is the furthest right in the chromosome, setting to None")
			self.right = None
		else: 
			self.right = gene_to_field[gene_i+1]
	
	
	def share_neighbouring_seqs(self):

		""" Determine which orientations the neighbouring seqs are in: head-to-head? head-to-tail? tail-to-head? tail-to-tail? """

		print("\tAttempting to share sequences between me and my leftward neighbour...")

		if self.left==None:
			# @TODO: make sure to deal with left-edge case
			pdb.set_trace()

		if self.right==None:
			# @TODO: make sure to deal with right-edge case
			pdb.set_trace()
		
		# Share seqs with left, Head-to-head (1/2 each), Head-to-tail (2/3 each)
		# @TODO: can reduce code by doing an IF strand directions are not the same then share 50:50 ...
		if ((self.left["c7_strand"]=="+") and (self.me["c7_strand"]=="-")) or ((self.left["c7_strand"]=="-") and (self.right["c7_strand"]=="+")):
			print("\t\tHead-to-Head (or Tail-to-Tail) case encountered: 5'===left===>3'......3'<===me===5'")	
			# @TODO: I assume the gff convention is to name the left-most chromosome pos as 0
			intergenic_seq_diff = self.me["c4_start"] - self.left["c5_end"]  # @TODO: code repetition, can factor out
			print("\t\t\tNo. intergenic BPs between me and left: "+str(intergenic_seq_diff))
			
			fraction_me = (1/2)
			print("\t\t\tFraction of BPs taken by me: "+str(fraction_me)) 
			
			me_share   = intergenic_seq_diff * fraction_me
			print("\t\t\tMy share of BPs: "+str(me_share)) 
			
			left_share = intergenic_seq_diff * (1-fraction_me)
			print("\t\t\tLeft share of BPs: "+str(left_share)) 

			# @LATEST:@TODO:add the shared intergenic positions to the self.me, or directly to file?
			
			# @TEST:@DONE: it ^, done: indeed correct orientation found and also correct sharing fractions

			self.me

		if ((self.left["c7_strand"]=="+") and (self.me["c7_strand"]=="+")):
			print("\t\tHead-to-Tail case encountered: 5'===left===>3'...intergenic...5'===me===>3'")

			intergenic_seq_diff = self.me["c4_start"] - self.left["c5_end"]
			fraction_me = (2/3) # @TODO: check with Kathrin: 2/3 to 5' or 2/3 to the 3'?
			me_share   = intergenic_seq_diff * fraction_me
			left_share = intergenic_seq_diff * (1-fraction_me)
			# @TEST: it^

		if ((self.left["c7_strand"]=="-") and (self.me["c7_strand"]=="-")):
			print("\t\tTail-to-Head case encountered: 3'<===left===5'...intergenic...3'<===me===5'")

			intergenic_seq_diff = self.me["c4_start"] - self.end["c5_end"]
			fraction_me = (1/3)
			me_share = intergeneic_seq_diff * fraction_me 
			left_share = intergenic_seq_diff * (1-fraction_me)
			# @TEST: it^
	
		# @LATEST
		# Share seqs with right, then take the share




# 3. We want to add two new fields to the gene_to_field dict, so that each gene has a "reg_start", and a "reg_end"

 

a = gene_and_neighbours( 2 )

a.share_neighbouring_seqs()


# }} 2 alternative 

