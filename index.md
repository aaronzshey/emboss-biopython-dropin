# Biopython as a drop-in replacement to EMBOSS

### Abstract

DNA sequencing remains one of the greatest innovations of the millenium, evolving from billion dollar behemoths such as the Human Genome Project to the modern-day Flongle, a 
device which fits in the palm of your hand and costs a mere $1,500.  Accompanying this rapid miniaturization of DNA sequencing hardware is the democratization of DNA 
sequencing software - where biologists used to need to purchase expensive software to analyze their sequences, they can now use myriad websites, apps, or programs to 
organize and understand their sequences.  However, software is never finished, dogged by the threat of a better solution just over the horizon.  

Since the early 2000s, the software package EMBOSS (European Molecular Biology Open Source Software) was the application *du jour*, used by labs all around the world. 
However, development of EMBOSS ceased in 2013, leaving a power vacuum soon filled by the package BioPython.  Migrating from EMBOSS to Biopython remains difficult and 
at times confusing.


### Introduction

In "A Python script to merge Sanger sequences," Chen et. al. present a simple and free Python program to merge multiple Sanger sequence results into a complete genome.
Running Chen et. al's program today returns multiple warnings, warning of outdateed EMBOSS and Biopython methods.  To ensure the continued function of this Python program, these outdated methods will be replaced by newer ones.


### Procedure

The offending block of code from Chen et. al's program occurs from lines 114 to 190:

```python
#function for the boundaries
"""
The function below parses the output of the EMBOSS needle alignment between the tandemly
arranged Sanger sequencing files to obtain the positions of the first 50 nucleotides match.
The EMBOSS needle shows the alignment result in segments of 50 nucleotides.
This boundary is used to determine which part of each Sanger sequence is kept for joining
to the final merged sequence. This operation is actually the same as that when we merge the
Sanger sequencing files manually.
"""
def needle_align_to_get_boundaries(f1,f2):
    alignment_a_left = 0
    alignment_b_left = 0
    output_file_name = f1.split('.')[0] + f2.split('.')[0] + ".needle"
    
    needle_cline = NeedleCommandline()
    needle_cline.asequence = f1
    needle_cline.bsequence = f2
    needle_cline.gapopen = 10
    needle_cline.gapextend= 0.5
    needle_cline.outfile = output_file_name
    print(needle_cline)
    stdout, stderr = needle_cline()
    print(stdout + stderr)

    #open the needle alignment output file and get boundaries
    file = open(output_file_name)
    file_lines = file.readlines()
    file.close()

    for line in file_lines:
        print(line, end="")

    alignment_a_squence_positions = []
    alignment_b_squence_positions = []

    file = open(output_file_name)

    new_line1 = file.readline()
    new_line2 = file.readline()

    while len(new_line2):
        line_a = new_line1
        line_b = new_line2
        new_line2 = new_line2.strip()
        
        if (50*'|' in new_line2):

            line_b = file.readline()
      
            alignment_a_squence_line_str = line_a.strip()
            alignment_b_squence_line_str = line_b.strip()
            print("The beginning of excellent alignment is shown below.\n")
            
            alignment_a_squence_line_str_split = alignment_a_squence_line_str.split()
            print(alignment_a_squence_line_str_split[0].ljust(5,' '),\
                  alignment_a_squence_line_str_split[1],\
                  alignment_a_squence_line_str_split[2].rjust(6,' '),\
                  sep="")
            alignment_b_squence_line_str_split = alignment_b_squence_line_str.split()
            print(alignment_b_squence_line_str_split[0].ljust(5,' '),\
                  alignment_b_squence_line_str_split[1],\
                  alignment_b_squence_line_str_split[2].rjust(6,' '),\
                  sep="")
     
            print("\n")

            alignment_a_left = int(alignment_a_squence_line_str.split()[0])
            alignment_b_left = int(alignment_b_squence_line_str.split()[0])
            break
        else:
            new_line1 = new_line2         #notice the skill here, we must go step by step through the lines
            new_line2 = file.readline()
    file.close()

    return alignment_a_left, alignment_b_left

#end of function
```

This code uses EMBOSS's `needle` command-line tool to align two sequence files, write the output to a file, and then read the file back to Python.  When the code
encounters a stretch of 50 vertical bars, or 50 perfect alignments, it copies the left-most coordinates (where the perfect alignments begin) and returns them for future use in the rest of the program.

To bring this code to the modern age, we must first use Biopython's `PairwiseAligner` function - a modern, more supported sequence aligner with support for a wider 
range of functions.  We must configure this pairwise aligner to act like `needle` by telling it to use the `EDNAFULL` or `NUC.4.4` substitution matrix and setting 
the gap open and extend values correctly.  Finally, we must find the leftmost coordinates of the beginnings of long chunks of alignments. 

This is one area where Biopython's Pairwise aligner is superior to EMBOSS' needle.  Instead of having to write our results to a file and then read it back, we can instead call the `.aligned` method to find the start and end coordinates of long chunks of alignments.  By finding the longest chunk of alignments, we can find the left-most positions of these alignments programatically.  

Here is the full code:

```python
def needle_align_to_get_boundaries(f1, f2):
    asequence = Path(f1).read_text().strip()
    bsequence = Path(f2).read_text().strip()
    
    # Define valid characters
    valid_chars = set("ACGTN")

    # Clean sequences
    asequence = ''.join([char for char in asequence if char in valid_chars])
    bsequence = ''.join([char for char in bsequence if char in valid_chars])

    # Initialize the aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10 
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = nuc44_matrix

    # Perform the alignment
    alignments = aligner.align(asequence, bsequence)

    '''
    print(alignment)
    GAACT
    ||--|
    GA--T

    alignment.aligned
    (((0, 2), (4, 5)), ((0, 2), (2, 3)))

    note how this tells us the fifth element of the target sequence (4, 5)
    is aligned to the third element of the query sequence (2, 3)

    We find the tuple with the greatest difference, representing the longest 
    chunk of perfect alignment.

    We also arbitrarily pick the first alignment, because it probably is the best one
    '''
    aligned = alignments[0].aligned

    # Convert alignments into diffs
    coord_diffs = list(map(lambda x: int(x[1] - x[0]), aligned[0]))
    # Find index of greatest diff in coords, indicating longest excellent alignment chunk
    max_diff_index = coord_diffs.index(max(coord_diffs))

    # Get left bounds from target (first aligned array) and query (second aligned array)
    alignment_a_left = aligned[0][max_diff_index][0]
    alignment_b_left = aligned[1][max_diff_index][0]
    
    return alignment_a_left, alignment_b_left
```

When pasting this function into Chen et. al's program, it functions exactly the same, in fewer lines of code.  

### Discussion
As the old adage goes, if it isn't broken, don't fix it.  While applicable to New Coke and Crystal Pepsi, it is less useful in terms of software.  Technical debt occurs when code becomes old, outdated, or unsupported - while old programs may still work fine, their continued function is not guaranteed, and upgrading their code will ensure their reliability.  

### Areas for Improvement
One glaring area for improvement is picking an alignment to use.  By default, Biopython's pairwise aligner returns the alignment with the highest score first - so we
pick the first alignment.  It is possible there is a better alignment somewhere in the many tens of thousands of alignments `PairwiseAligner` finds. 
