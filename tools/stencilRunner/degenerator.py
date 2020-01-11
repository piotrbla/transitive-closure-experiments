with open('decompiler.bat', 'w') as fdeco, open('testrun.bat', 'w') as frun:
	for i in range(1, 400, 2):
		s= str(i+1)
		out_filename = 'nont_' + s + '.cpp' 
		copy_file = False
		with open('nontiled.cpp') as f, open(out_filename, 'w') as fout, open('cDYN0_' + s + '.out' ) as fin:
			for line in f:
				if line.startswith('#insertcode'):
					for linein in fin:
						if linein.startswith('for (int') or copy_file: 
							if linein.startswith('End'):
								break
							fout.write(linein)
							copy_file = True
					print(s)
				else:
					fout.write(line)
		fdeco.write('cl /GS /GL /analyze- /W3 /Gy /Zc:wchar_t /Zi /Gm- /O2 /sdl /Fd /Zc:inline /fp:precise /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_UNICODE" /D "UNICODE" /errorReport:prompt /WX- /Zc:forScope /Gd /Oy- /Oi /MD /FC /nologo  /diagnostics:column ' + out_filename + '\n')
		frun.write('nont_' + s + '.exe' '\n')
