for i in range(1, 70, 2):
	s= str(i+1)
	print("printf \"26\\n30\\n" + s + "\\n\" > tile.sizes")
	print("./polycc --tile --parallel NusValidation.cpp")
	print("gcc -o nus" + s + " -O2 NusValidation.cpp.pluto.c -fopenmp -lm")
	print("cp nus" + s + " exp")
	print("./exp/nus" + s)
	

