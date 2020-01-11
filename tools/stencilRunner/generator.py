for i in range(1, 400, 2):
	s= str(i+1)
	print("..\stencil.exe filename=./computeDYN0.c sizes=3 size1=16 size2=" + s + " size3=64 method=1 > cDYN0_" + s + ".out 2> cDYN0_" + s + ".dbg")