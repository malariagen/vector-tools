import sys
import numpy as np

sel=[]
list_files = sys.argv[1].replace(" ","")[1:-1].split(',')
for i in list_files:
	if sys.argv[2] in i:
		sel.append(i)
print(sel)
