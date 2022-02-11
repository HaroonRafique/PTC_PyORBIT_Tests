#
# Python script to split the downloaded file from the LayoutDB
# and generate the various files:
#  - psb.ele  lines between ELEMENT - END_ELEMENT in the input file
#  - psb.seq  lines between SEQUENCE - END_SEQUENCE in the input file
#  - psb.str  lines between STRENGTH - END_STRENGTH in the input file
#
# (C) IE - 31.03.2017
#

import argparse

parser = argparse.ArgumentParser(description='Split MAD-X generated file from Layout DB')
parser.add_argument('-i','--inputfile', dest='inputfile', required="True")
parser.add_argument('-o','--outputfile', dest='outputfile', required="True")

args = parser.parse_args()
inputfile = args.inputfile
outputfile = args.outputfile

print 'Input file is "', inputfile
print 'Output file root is "', outputfile

head = []
ihead = 0
iele = 0
iseq = 0
istr = 0

input_file = open(inputfile, 'r')
count_lines = 0
count_head = 0
count_ele = 0
count_seq = 0
count_str = 0

out_ele = open(outputfile+'.ele','w')
out_seq = open(outputfile+'.seq','w')
out_str = open(outputfile+'.str','w')

for line in input_file:
    #  print line
    count_lines += 1
    if "++HEAD++" in line:
       ihead = 1
       print line
    if "++END_HEAD++" in line:
       ihead = 0
       out_ele.write("\n\n")
       out_seq.write("\n\n")
       out_str.write("\n\n")
    if ihead == 1:
       count_head += 1
       if count_head > 1 :
           out_ele.write(line)
           out_seq.write(line)
           out_str.write(line)
           head.append(line)

    if "++ELEMENTS++" in line:
       iele = 1
       print line
    if "++END_ELEMENTS++" in line:
       iele = 0
    if iele == 1:
       count_ele += 1
       if count_ele > 1:
           out_ele.write(line)
    if "++SEQUENCE++" in line:
       iseq = 1
       print line
    if "++END_SEQUENCE++" in line:
       iseq = 0
    if iseq == 1:
       count_seq += 1
       if count_seq > 1:
           out_seq.write(line)

    if "++STRENGTH++" in line:
       istr = 1
       print line
    if "++END_STRENGTH++" in line:
       istr = 0
    if istr == 1:
       count_str += 1
       if count_str > 1:
           out_str.write(line)

# print 'header:', head
print 'number of lines           :', count_lines
print 'number of header lines    :', count_head
print 'number of element lines   :', count_ele
print 'number of sequence lines  :', count_seq
print 'number of srtrength lines :', count_str
print ' total = ', count_head + count_ele + count_seq + count_str

out_ele.close()
out_seq.close()
out_str.close()
