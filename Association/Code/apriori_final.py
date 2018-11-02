import os
import csv
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA 
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import itertools

s = int(input("Enter the support threshold: "))
min_conf = int(input("Enter the confidence threshold: "))


def cal_support_count(srows, threshld):
	return (threshld * srows)/100

def template1_func(query,rules_set):
	part1=query[0]
	part2=query[1]
	part3=(query[2]).replace('[','').replace(']','').split(',')
	pruned_ruleset=[]
	for rule in rules_set:
		rulelist=rule.split('->')
		head_list=rulelist[0].split(',')
		body_list=rulelist[1].split(',')
		mainset=[]
		if(part1=='RULE' and part2=='ANY'):
			mainset = set(head_list).union(set(body_list))
			if(len(mainset & set(part3))>0):
				pruned_ruleset.append(rule)
			# print(pruned_ruleset)

		elif(part1=='RULE' and part2=='NONE'):
			mainset = set(head_list).union(set(body_list))
			if(len(mainset & set(part3))==0):
				pruned_ruleset.append(rule)
			
		elif(part1=='RULE' and (part2 is not 'ANY' and  part2 is not 'NONE')):
			mainset = set(head_list).union(set(body_list))
			if(len(mainset & set(part3))==int(part2)):
				pruned_ruleset.append(rule)
			
		elif(part1=='HEAD' and part2=='ANY'):
			mainset=set(head_list)
			if(len(mainset & set(part3)) >0):
				pruned_ruleset.append(rule)
			
		elif(part1=='HEAD' and part2=='NONE'):
			mainset=set(head_list)
			if(len(mainset & set(part3)) ==0):
				pruned_ruleset.append(rule)

		elif(part1=='HEAD' and (part2 is not 'ANY' and  part2 is not 'NONE')):
			mainset=set(head_list)
			#print(mainset)
			if(len(mainset & set(part3)) ==int(part2)):
				pruned_ruleset.append(rule)  
			 

		elif(part1=='BODY' and part2=='ANY'):
			mainset=set(body_list)
			if(len(mainset & set(part3)) >0):
				pruned_ruleset.append(rule)
		elif(part1=='BODY' and part2=='NONE'):
			mainset=set(body_list)
			if(len(mainset & set(part3)) ==0):
				pruned_ruleset.append(rule)

		elif(part1=='BODY' and (part2 is not 'ANY' and  part2 is not 'NONE')):
			mainset=set(body_list)
			if(len(mainset & set(part3)) ==int(part2)):
				pruned_ruleset.append(rule) 
	# print(pruned_ruleset)
	return pruned_ruleset 

def template2_func(query,rules_set):
	part1=query[0]
	part2=query[1]
	#print(part1,part2)
	pruned_ruleset=[]

	for rule in rules_set:
		rulelist=rule.split('->')
		head_list=rulelist[0].split(',')
		body_list=rulelist[1].split(',')
		mainset=[]
		if(part1=='RULE'):
			mainset=set(head_list).union(set(body_list))
			if(len(mainset)==int(part2)):
				pruned_ruleset.append(rule)
			#print(pruned_ruleset)
		elif(part1=='HEAD'):
			mainset=set(head_list)
			if(len(mainset)==int(part2)):
				pruned_ruleset.append(rule)
			#print(pruned_ruleset)
		elif(part1=='BODY'):
			mainset=set(body_list)
			if(len(mainset)==int(part2)):
				pruned_ruleset.append(rule)
	#print(pruned_ruleset)
	return pruned_ruleset


with open('associationruletestdata.txt', 'r') as in_file:
	stripped = (line.strip() for line in in_file)
	lines = (line.split(",") for line in stripped if line)
	with open('asstestdata.csv', 'w') as out_file:
		writer = csv.writer(out_file)
		writer.writerows(lines)

def get_itemset_rules(itemset_list, body_max_length, numerator, confidence_threshold_percentage, gene_matrix, rules_set):
	itemset_list_set = set(itemset_list)

	for body_length in range(body_max_length, 0, -1):
		for body_list in list(itertools.combinations(itemset_list, body_length)):
			body_set = set(body_list)
			head_list = list(body_set ^ itemset_list_set)
			rule = ",".join(body_list) + "->" + ",".join(head_list)
			if rule not in rules_set:
				denominator = 0
				for gene_list in gene_matrix:
					gene_list_set = set(gene_list)
					if body_set < gene_list_set:
						denominator += 1
				confidence_percentage = (numerator / denominator) * 100
				if (confidence_percentage >= confidence_threshold_percentage):
					rules_set.add(rule)
	return rules_set


gdata = []


apriori = pd.read_csv("asstestdata.csv",sep='\t',header = None)

rows = apriori.shape[0]
columns	= apriori.shape[1]

disease = apriori.shape[1]

for i in range(disease-1):
	gdata.append('G'+ str(i))

gdata.append('disease')

apriori.columns = gdata

for i in range(0,apriori.shape[1]-1):
	for j in range(0,apriori.shape[0]):
		apriori.iloc[j,i] = 'G' + str(i) + str(apriori.iloc[j,i])

superlist = []


for i in range(0,rows):
	superlist.append( [apriori.values[i,j] for j in range(columns-1) ] )

srows = len(superlist)
scols = len(superlist[0])


frq_items1 = Counter()

for i in range(0,scols):
	for j in range(0,srows):
		frq_items1[superlist[j][i]] += 1



sup_count = cal_support_count(srows, s)
#print(sup_count)
validated_frq_items1 = Counter()

############################ VALID FREQUENT ITEM SENT PRINT ########################################33
confid_map = Counter()

confid_map = frq_items1

for item in frq_items1:

	if frq_items1[item] >= sup_count:
		validated_frq_items1[item] = frq_items1[item]
# print(len(validated_frq_items1))

k = 2
k_frq_list = []
k_frq_list.append((validated_frq_items1))

valid_freq_list = []
valid_freq_list.append(validated_frq_items1)

while (1):
	k_touple = (list(itertools.combinations(valid_freq_list[k-2],k)))	

	for values in k_touple:
		touple = []
		for i in range(0,len(values)):
			touple.append(int(str(values[i]).split('G')[1].split('D')[0].split('U')[0]))

		#remove G1_UP and G!Down cases
		if len(set(touple)) == 1:
			k_touple.remove(values)


	valid_set2 = Counter()

	for item in k_touple:
		all_items = []
		for i in range(0,len(item)):
			all_items.append(int(str(item[i]).split('G')[1].split('D')[0].split('U')[0]))

		c1 = 0

		for i,row in enumerate(superlist):
			values = 1
			for j,col in enumerate(all_items):
				values = values & (superlist[i][col] == item[j])

				if values == 0:
					break
			c1+=values

		valid_set2[item] = c1

	validated_frq_items2 = Counter()
	#print(len(valid_set2))

	
	for key,v in valid_set2.items():
		#print(ord(key[0]),ord(key[0]))
		#int(str(item[i]).split('G')[1].split('D')[0].split('U')[0])
		
		#key.sort()
		confid_map[tuple(sorted(key))] = v
	#print(confid_map)

	#minimum support items retain	
	cnt = 0
	for item in valid_set2:

		if valid_set2[item] >= sup_count:
			validated_frq_items2[item] = valid_set2[item]
			cnt+=1

	k_frq_list.append((validated_frq_items2) )
	
	valid_list = []
	for key in validated_frq_items2.keys():
		for t,val in enumerate(key):
			valid_list.append(key[t])

	valid_items = set(valid_list)

	valid_freq_list.append(valid_items)
#############################PRINT FREQUENT ITEM SET AND LENGTH #################################################
	# print(len(validated_frq_items2),cnt)

	k+=1
	if cnt == 0:
		break

#######################		RULE GENERATION########################################################################
#print(k_frq_list[2].keys())
dummy_set = set()
rules_set = set()
for i in range(1,len(k_frq_list)):
	
	for touple in k_frq_list[i].keys():
		#set_touple = set(touple)
		maxlen = len(touple) - 1

		for k in range(maxlen,0,-1):
			
			for entry in list(itertools.combinations(touple,k)):

				#print(entry)
				#set_entry = set(entry)
				#print(set_entry)
				head_list = list(set(entry) ^ set(touple))

				set_entry_tuple = []
				head_list_tuple = []
				rule = ",".join(entry) + "->" + ",".join(head_list)

				touple = tuple(sorted(touple))
				if len(entry) == 1:
					entry = entry[0]
				else:
					entry = tuple(sorted(entry))
				#print(entry)
				dummy_set.add(rule)
				# for s in entry:
#				print(entry)
				confidence = (confid_map[touple] / confid_map[entry]) * 100

				#confidence = k_frq_list[k][item] / confid_map[item[0]] * 100
				if confidence >= min_conf:
					rules_set.add(rule)



#print(len(dummy_set))
# print("rules",rules_set)
# print("rules_set",len(rules_set))

#print(rules_set)
#sys.exit()
'''
============================================ Template1=========================

("RULE", "ANY", ['G71_UP','G81Down'])

================================================================================


============================================ Template2=========================

("RULE", 3)

================================================================================




============================================ Template3=========================

1and1|RULE:ANY:[G71Up,G81Down];RULE:ANY:[G71Up,G81Down]

1and2|RULE:ANY:[G71Up,G81Down];BODY:2

2and2|RULE:3;BODY:2

================================================================================


'''



template = input('Enter template type: ')
if(template=='1'):

	query=input('Enter template-1 query ').replace('_','').replace('P','p').replace(' ','').replace('"','').replace("'",'').replace('(','').replace(')','').replace(',',':').split('[')
	
	#print(query)
	first=query[0]
	#print(first)
	second=query[1].replace(':',",")
	second="["+second
	#print(second)
	query=first+second
	#print(query)
	query=query.split(":")
	#print(query)



	query=tuple(query)
	pruned_rule_set=template1_func(query,rules_set)
	print(pruned_rule_set)
	print(len(pruned_rule_set))
#print(pruned_ruleset)



if(template=='2'):
	query=input('Enter template-2 query:').replace('"','').replace("'",'').replace('(','').replace(')','').replace(' ','').replace(",",":").split(':')
	
	query=tuple(query)
	pruned_rule_set=template2_func(query,rules_set)
	# print(pruned_rule_set)
	print(pruned_rule_set)
	print(len(pruned_rule_set))
# 1and1|RULE:ANY:[G71Up,G81Down];RULE:ANY:[G71Up,G81Down]
if(template=='3'):
	query=input('Enter template-3 query:').split('|')
	query=tuple(query)
	format1=query[0]
	format2=query[1]
	templates=format2.split(';')
	templates[0]=tuple(templates[0].split(":"))
	templates[1]=tuple(templates[1].split(":"))

	#print(templates[0],templates[1])




	if(format1[0]=='1' and format1[-1]=='1' and 'or' in format1):
		pruned_rule_set=(set(template1_func(templates[0],rules_set)).union(set(template1_func(templates[1],rules_set))))
		print(pruned_rule_set)
	if(format1[0]=='1' and format1[-1]=='1' and 'and' in format1):
		pruned_rule_set=(set(template1_func(templates[0],rules_set)).intersection(set(template1_func(templates[1],rules_set))))
	
	if(format1[0]=='1' and format1[-1]=='2' and 'or' in format1):
		pruned_rule_set=(set(template1_func(templates[0],rules_set)).union(set(template2_func(templates[1],rules_set))))

	if(format1[0]=='1' and format1[-1]=='2' and 'and' in format1):
		pruned_rule_set=(set(template1_func(templates[0],rules_set)).intersection(set(template2_func(templates[1],rules_set))))

	if(format1[0]=='2' and format1[-1]=='2' and 'or' in format1):
		pruned_rule_set=(set(template2_func(templates[0],rules_set)).union(set(template2_func(templates[1],rules_set))))

	if(format1[0]=='2' and format1[-1]=='2' and 'and' in format1):
		pruned_rule_set=(set(template2_func(templates[0],rules_set)).intersection(set(template2_func(templates[1],rules_set))))	

	print(pruned_rule_set)
	print(len(pruned_rule_set))
    







		



