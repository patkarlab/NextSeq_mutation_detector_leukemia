import sys
import csv
#combine = sys.argv[1]
#somatic = sys.argv[2]
#outfile = sys.argv[3]

mergeList = []
finalList = []
#with open(sys.argv[3], 'w') as outfile:
with open(sys.argv[1], 'r') as combine:
	for cline in combine:
		#print(cline.rstrip())
		clist = [cline.rstrip().split(',')]
		#print(clist)
		c = 0
		with open(sys.argv[2], 'r') as somatic:
			for sline in somatic:
				slist = [sline.rstrip().split(',')]
				check_slist = 0
				#print(slist)
				if clist[0][0] == slist[0][0] and clist[0][1] == slist[0][1]:
					mergeList.append(slist[0][:2])
					slist[0][5] = slist[0][5] + clist[0][5]
					num = int(slist[0][8]) + int(clist[0][8])
					slist[0][8] = num
					finalList.append(slist[0])
					#outfile.write(str(slist))
					c = 1
					#print("match")	
				num = 0
		if c == 0:
			#outfile.write(str(clist))
			finalList.append(clist[0])
			mergeList.append(clist[0][:2])
#print(mergeList)
with open(sys.argv[2], 'r') as somatic:
	for sline in somatic:
		slist = [sline.rstrip().split(',')]
		if slist[0][:2] not in mergeList:
			#outfile.write(str(sline))
			finalList.append(slist[0])

with open(sys.argv[3], 'w', newline='') as outfile:
	writer = csv.writer(outfile)
	writer.writerows(finalList)

#with open(sys.argv[3], 'w') as outfile:
#	for item in finalList:
#		#outfile.write("%s\n" % item)
#		outfile.write(item + "\n")
outfile.close()		
