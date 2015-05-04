#!/usr/bin/env python
# encoding: utf-8
"""
sam2bed.py

Created by Aaron Quinlan on 2009-08-27.
Copyright (c) 2009 Aaron R. Quinlan. All rights reserved.
"""

import sys
import getopt
import re

help_message = '''
	USAGE: sam2bed -s <sam>
'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def processSAM(file):
	"""
		Load a SAM file and convert each line to BED format.
	"""		
	f = open(file,'r')
	for line in f.readlines():
		samLine = splitLine(line.strip())
		makeBED(samLine)
	f.close()	
	
					
def makeBED(samFields):
	
	samFlag = int(samFields[1])
	
	# Only create a BED entry if the read was aligned
	if (not (samFlag & 0x0004)):
		
		chrom = samFields[2]
		start = str(int(samFields[3])-1)
		end = str(int(samFields[3]) + len(samFields[9]) - 1)
		name = samFields[0]	
		strand = getStrand(samFlag)

		# Let's use the edit distance as the BED score.
		# Clearly alternatives exist.
		editPattern = re.compile('NM\:i\:(\d+)')
		editDistance = editPattern.findall(samFields[12])

		# Write out the BED entry
		print chrom + "\t" + start + "\t" + end + "\t" + name + "\t" + editDistance[0] + "\t" + strand
		
		
def splitLine(line, delim="\t"):
	splitline = line.split(delim)
	return splitline		


def getStrand(samFlag):
	strand = "+"
	if (samFlag & (0x10)):	# minus strand if true.
		strand = "-"		
	return strand
