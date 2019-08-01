import sys
import pyfasta

fa = pyfasta.Fasta(sys.argv[1], key_fn=lambda key: key.split()[0]) 

d = {}
for k in fa.keys():
	d.update({k.split()[0]: fa[k]})

print(len(d[sys.argv[2]]))
