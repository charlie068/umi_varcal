
import os
import sys
import glob

thisProgramName  = sys.argv[0]

try:
	RAMScriptName   = sys.argv[1]
except:
	RAMScriptName = "RAM.py"

try:
	setupScriptName  = sys.argv[2]
except:
	setupScriptName = "setup.py"


scripts = glob.glob('*.py')
scripts.remove(thisProgramName)
scripts.remove(RAMScriptName)
scripts.remove(setupScriptName)
scripts.remove("TreatReads.py")






pyxFile = open("functions.pyx", "w")
imports = ""
code = ""


for script in scripts:
	scripFile = open(script)
	for line in scripFile:
		if "import" in line:
			line = line.replace("\t", "")
			if line not in imports and "*" not in line and "#" not in line:
				imports += line.replace("\t", "")
		else:
			code += line

	code += "\n\n"



pyxFile.write(imports+"\n\n"+code+"\n\n")

print("\nDone\n")