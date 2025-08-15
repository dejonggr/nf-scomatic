#!/usr/bin/env python3
import pysam
import sys
import os

#run with the bam name as the first positional argument
#(and have the bam name be the sample name)
bampath = sys.argv[1]
sample = os.path.splitext(bampath)[0]

#open up the original bam, and the new one which will have CB flagged
bamfile = pysam.AlignmentFile(bampath, "rb")
tweaked = pysam.AlignmentFile(sample+"-flagged.bam", "wb", template=bamfile)
for read in bamfile.fetch():
    #do we have a CB tag? if so, stick the sample on as a prefix. write out to flagged
    if read.has_tag('CB'): 
        #hooray, weird starsolo bams with - CBs sometimes! yay!
        if read.get_tag('CB') == '-':
            #wipe useless tags to avoid confusing
            read.set_tag('CB', None)
            read.set_tag('UB', None)
        else:
            #normal read, tag properly
            read.set_tag('CB',sample+'-'+read.get_tag('CB'))
    tweaked.write(read)

tweaked.close()
bamfile.close()
