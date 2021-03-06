with open('testrun.sh', 'w') as frun, open('copier.bat', 'w') as fcopier:
	fcopier.write('copy testrun.sh dist\\testrun.sh\n')
	for i in range(11, 20, 2):
		s= str(i+1)
		out_filename = 'tiled_' + s + '.cc' 
		out_exe_filename = 'tiled_' + s 
		copy_file = False
		with open('tiled.cc') as f, open(out_filename, 'w') as fout, open('cBST_' + s + '.out' ) as fin:
			for line in f:
				if line.startswith('#insertcode'):
					for linein in fin:
						if linein == '{\n' or linein.startswith('for (int') or copy_file: 
							if linein.startswith('End'):
								break
							fout.write(linein)
							copy_file = True
					print(s)
				else:
					fout.write(line)
		frun.write('.\\' + out_exe_filename + '\n')
		fcopier.write('copy ' + out_filename + ' dist\\' + out_filename +  '\n')
		