#name: Duy Truong
#class: Bioinformatic

import random
from collections import OrderedDict


f = open('CLIB215_AEWP01000000_cds.fsa', 'r')
o = open('CLIB215_AEWM01000000_cds_del.fsa', 'w')
lines = f.readlines();
stri = ""
for line in lines:
	if line[0] != '>':
		line = line.replace('\n', '')
		stri += line
o.write(stri)
o.close()
f.close()

length = 150;
k_mer = 5;
f = open('CLIB215_AEWM01000000_cds_del_A.fsa', 'r')
n = len(f.read())
#get the DNA sequences and get random fragment
stop= n-length
times = n/length
f.seek(0,0)
newline=f.read()
newline = newline.replace("\n","")
original = newline
array = [0 for i in range(0,times)]
for i in range (0, times):
	jump = random.randint(0,stop)
	array[i] = newline[jump:jump+length]
f.close()
#keep track the over lap
list(OrderedDict.fromkeys(array))
o = open("overlap2.fna","w")
overlap = [[0 for m in range(0,times)] for m in range(0,times)]

for m  in xrange(0,times):
	for n  in xrange(0,times):
		if(m==n):
			overlap[m][m]= 0
		check = array[m][length-k_mer-1: length-1]
		check1= array[n][0:k_mer]
		if(check1 == check):
			overlap[m][n]= len(check1)
			
		o.write(str(m)+", " +str(n)+", " +str(overlap[m][n])+"\n")
o.close()



newDNA = [0 for m in range(0,times)]
for m in range(0,times):
	newDNA[m] = array[m]

biggest_overlap= [0 for m in range(0,times)]
index_biggest =  [0 for m in range(0,times)]


for x in xrange(0,times):
	index = 0
	biggest_overlap[x] = max(overlap[x])
	while (overlap[x][index] != max(overlap[x])):
		index+=1
	index_biggest[x] = index
#merge the DNA gradment to get new DNA
temp = biggest_overlap
temp1 = index_biggest
for x in xrange(0,times):
	m=x
	while (temp [index_biggest[m]] != 0):
		newDNA[x] += array[index_biggest[m]][biggest_overlap[m]+1: length-1]
		temp [index_biggest[m]] =0
		m = index_biggest[m]
print newDNA

o = open("overlap_original1.fna","w")
#compare with the original chromosomes				
for x in xrange(0, times):
	original_overlap = " "
	for y in xrange (0, len(newDNA[x])):
		n = len(newDNA[x])
		for z in xrange (0, len(original)- n):
			l = len (newDNA[x][y:len(newDNA[x])] )
			if( newDNA[x][y:len(newDNA[x])] == original[z:z+l] ):
				temp = newDNA[x][y:len(newDNA[x])]
			if(len(temp) > len (original_overlap)):
				original_overlap = temp
	o.write("newDNA number: "+str(x)+"\n" + "length: " + str(len(original_overlap)) + "\n" +"overlap with original: " + str(original_overlap) + "\n")
	print original_overlap
o.close()
