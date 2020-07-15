# Copyright (C) 2016 VU University Medical Center Amsterdam
# Author: Roy Straver (r.straver@vumc.nl)
#
# This file is part of SANEFALCON
# SANEFALCON is distributed under the following license:
# Attribution-NonCommercial-ShareAlike, CC BY-NC-SA (https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
# This license is governed by Dutch law and this license is subject to the exclusive jurisdiction of the courts of the Netherlands.



import sys
from collections import OrderedDict

reference_file_name = sys.argv[1]
profile_file_name = sys.argv[2]

nuclData=dict()
with open(reference_file_name, 'r') as referenceFile:
	correlations = [float(x) for x in referenceFile.readline().split(" ")]
	scalars = [float(x) for x in referenceFile.readline().split(" ")]
#print len(correlations),scalars
	
with open(profile_file_name, 'r') as profileFile:
	# In these readlines -1 to drop the 'end of line' comma from a previous awk
	my_lines = profileFile.readlines()
profileFile.close()
profile = [float(x) for x in my_lines[0].split(",") if x != '\n']
profile.reverse()

profile.extend([float(x) for x in my_lines[1].split(",") if x != '\n'])

no_dup_profile = list(OrderedDict.fromkeys(profile))
try:
	assert len(correlations) == len(no_dup_profile)
except AssertionError:
	print('Data not aligned on ', reference_file_name)
	exit(1)

totalReads = sum(no_dup_profile)
normProfile = [x/totalReads for x in no_dup_profile]
	
summed = 0
for i, val in enumerate(normProfile):
	summed += val*correlations[i]
fetalFraction = scalars[0]*summed+scalars[1]

if len(correlations) == len(normProfile):
	print("Fetal Fraction:\t", fetalFraction)
else:
	print("ERROR: correlation and sample profiles are not aligned")
print("Nucleosome Profile:", "\t".join([str(x) for x in normProfile]))
