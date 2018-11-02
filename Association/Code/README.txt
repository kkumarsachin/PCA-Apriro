Run code as follows

python3.apriori_final.py

#Take input from command-line

Enter the support threshold: 50
Enter the confidence threshold: 70

#To print frequent item sets:

print(validated_frq_items2)

#To run association analysis based on template queries:

#Example Template 1
Enter template type: 1 
Enter template-1 query:("RULE", "ANY", ['G71_UP','G81Down'])

#Example Template 2
Enter template type: 2
Enter template-2 query: ("RULE", 3)


#Example Template 3
Enter template type: 3 
Enter template-3 query:1and2|BODY:ANY:[G1Down];HEAD:2
