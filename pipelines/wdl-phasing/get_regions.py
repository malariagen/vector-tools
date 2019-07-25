import sys

s=0
e=int(sys.argv[3])
while s <= int(sys.argv[1]):
	print(str(s) + "	" + str(e))
	s = e - int(sys.argv[2])
	e = s + int(sys.argv[3])

