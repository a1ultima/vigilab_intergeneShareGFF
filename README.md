
# VigiLab's intergeneShareGFF tool (for Kathrin)

---

## Description

Takes a GFF file, and extends the "START" and "END" fields of each gene row, to include intergenic regions

	data in: 
		<some GFF file>

	data out:
		<the same GFF file, with new START and END fields, to incorporate the intergenic regions>

**Here's a visual:**

We want to share the intergenic regions (i.e. "...") between, say, transcripts in a GFF file. Then print out a new GFF file with modified "start" and "end" fields to aknowledge these sharings. The following sketch illustrates all possible scenarios we need to account for, and the `|` shows the precise location along the intergenic regions, where nucleotide basepairs (BPs, shown as `.`) are shared between tandem gene/transcript/exon/cds annotations of the input GFF file.

| Case  | Diagram | Strands | Share  | Status  |
|---|---|---|---|---|
| Tail-to-Tail  | 5'(head)===>3'(tail)...\|...3'(tail)<===5'(head)  | (+,-)  | 1/2:1/2  | Tested  |
| Head-to-Head  | 3'(tail)<===5'(head)...\|...5'(head)===>3'(tail)  | (+,-)  | 1/2:1/2  | Tested  |
| Tail-to-Head  | 5'(head)===>3'(tail)..\|....5'(head)===>3'(tail)  | (+,+)  | 1/3:2/3  | Testing:@TODO  |
| Head-to-Tail  | 3'(tail)<===5'(head)....\|..3'(tail)<===5'(head)  |  (-,-) | 2/3:1/3  |  Developing:@TODO |

In the input GFF file, the above data is represented as:

`field: \/start\/end`

`row 1: ===1===>`

`row 2: <==2==`

`row 3: <==3==`

`row 4: ==4==>`

`row 5: ===5===>`

`row 6: <===6===`

Specifically, we want to distribute intergenic regions between tandem gene/transcript/exon/cds (@features), according to the following ratios: 1/2:1/2 and 1/3:2/3, inferred using the tandem gene's pairwise "strand" field configurations (see table above). And create a new GFF with modified "start" and "end" fields, which account for the, now shared, BPs (`.`) as shown above:

`row 1: ===1===>...`

`row 2: ...<==2==....`

`row 3: ..<==3==...`

`row 4: ...==4==>....`

`row 5: ..===5===>...`

`row 6: ...<===6===`

---

# REQUIREMENTS
 - Python3.5m
