import sys
import pyfasta

fa = pyfasta.Fasta(sys.argv[1]) 

print(len(fa[sys.argv[2]]))
