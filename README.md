
#DESCRIPTION

Takes a GFF file, and extends the "START" and "END" fields of each gene row, to include intergenic regions


	data in: 

		<some GFF file>

	data out:

		<the same GFF file, with new START and END fields, to incorporate the intergenic regions>
*Here's a visual:

We want to share the intergenic regions (i.e. "...") between, say, transcripts in a GFF file. Then print out a new GFF file with modified "start" and "end" fields to aknowledge these sharings. 

`===1====>...|...<==2==....|..<==3==...|...==4==>....|..====5====>...|...<===6===`

In the input GFF file, the above data is represented as:

`field: \/start\/end`

`row 1: ===1===>`

`row 2: <==2==`

`row 3: <==3==`

`row 4: ==4==>`

`row 5: ===5===>`

`row 6: <===6===`

But we want to aknowledge the shared intergeneic regions, so that intergenic regions are shared at a 1/2:1/2 and 1/3:2/3 ratio depending on the pairwise "strand" fields of tandem neighbouring transcripts, as such:

`row 1: ===1===>...`

`row 2: ...<==2==....`

`row 3: ..<==3==...`

`row 4: ...==4==>....`

`row 5: ..===5===>...`

`row 6: ...<===6===`

And create a new GFF with modified "start" and "end" fields, which account for 
