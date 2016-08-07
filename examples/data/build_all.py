
import os
import sys

print(sys.argv)

pdbdir = os.path.dirname(sys.argv[0])

print(pdbdir)

pdb_files = [f for f in os.listdir(pdbdir) if f[-4:] == '.pdb']

for gs in ["1.0", "0.5", "0.25"]
  for f in pdb_files:
    ex = "%s %s --grid-spacing %s" % (sys.argv[1], pdbdir + '/' + f, gs)
    print(ex)
    for line in os.popen(ex):
      print("!" + line)





