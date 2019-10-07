import sys
import numpy as np

sel=""
for f in sys.argv[2:]:
	if sys.argv[1] in f:
		sel += f + "	"
print(sel[:-1])
