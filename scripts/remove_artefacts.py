import sys
import csv

## arg1 => concat combined + somatic csv file
## arg2 => artefact csv
## arg3 => final csv without artefact outfile
## arg4 => artefact variants outfile

blackList = [] #artefact list
finalList = [] #without artefact final list

with open(sys.argv[1], 'r') as variantFile:
	for variant in variantFile:
		#print(cline.rstrip())
		list1 = [variant.rstrip().split(',')]
		#print(clist)
		counter = 0
		with open(sys.argv[2], 'r') as artefactFile:
			for artefact in artefactFile:
				list2 = [artefact.rstrip().split(',')]
				if list1[0][0] == list2[0][0] and list1[0][1] == list2[0][1] and list1[0][4] == list2[0][3]:
					blackList.append(list1[0])
					counter = 1
					#print("match")	
		if counter == 0:
			finalList.append(list1[0])

# writing lists to csv files
with open(sys.argv[3], 'w', newline='') as outfile_final:
	writer = csv.writer(outfile_final)
	writer.writerows(finalList)

with open(sys.argv[4], 'w', newline='') as outfile_artefacts:
	writer = csv.writer(outfile_artefacts)
	writer.writerows(blackList)

outfile_final.close()
outfile_artefacts.close()
