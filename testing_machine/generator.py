for i in range(11, 20, 2):
	s= str(i+1)
	print("..\dapt.exe filename=./nontiled.c sizes=3 size1=32 size2=16 size3=" + s + " method=1 > cBST_" + s + ".out 2> cBST_" + s + ".dbg")