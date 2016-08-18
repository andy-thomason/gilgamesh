
# build molecules in examples.
# takes one parameter: the path to the "molecules" executable.

import os
import sys

print(sys.argv)

pdbdir = os.path.dirname(sys.argv[0])

print(pdbdir)

pdb_files = [f for f in os.listdir(pdbdir) if f[-4:] == '.pdb']

if False:
  # list the chains
  for f in pdb_files:
    ex = "%s %s --list-chains" % (sys.argv[1], pdbdir + '/' + f)
    print(ex)
    for line in os.popen(ex):
      print(line)

chains = [
  ("2PTC.pdb", "E", "I"),
  ("1OYV.pdb", "AB", "I"),
  ("1OHZ.pdb", "A", "B"),
  ("4KC3.pdb", "A", "B"),
  ("1GRN.pdb", "A", "B"),
  ("1FSS.pdb", "A", "B"),
  ("1C4Q.pdb", "ABDE", "C"),
  ("1J1E.pdb", "ABCDE", "F"),
  ("1EMV.pdb", "A", "B"),
  #("4FOZ.pdb", "A", ""),
  ("1N8Z.pdb", "AB", "C"),
  ("2ZW3.pdb", "ABCDE", "F"),
  ("2VIR.pdb", "AB", "C"),
  ("1BRS.pdb", "ABCDE", "F"),
  ("4GRG.pdb", "ABC", "D"),
]

# thumbnail meshes
for f, ch1, ch2 in chains:
  ex = "%s %s se --chains %s%s --lod 0" % (sys.argv[1], pdbdir + '/' + f, ch1, ch2)
  print(ex)
  for line in os.popen(ex):
    print(line)

# LOD 1 solvent excluded models
for f, ch1, ch2 in chains:
  ex = "%s %s se --chains %s --lod 1" % (sys.argv[1], pdbdir + '/' + f, ch1)
  print(ex)
  for line in os.popen(ex):
    print(line)

  ex = "%s %s se --chains %s --lod 1" % (sys.argv[1], pdbdir + '/' + f, ch2)
  print(ex)
  for line in os.popen(ex):
    print(line)

# ball and stick models
for f, ch1, ch2 in chains:
  ex = "%s %s bs --chains %s --lod 1" % (sys.argv[1], pdbdir + '/' + f, ch1)
  print(ex)
  for line in os.popen(ex):
    print(line)

  ex = "%s %s bs --chains %s --lod 1" % (sys.argv[1], pdbdir + '/' + f, ch2)
  print(ex)
  for line in os.popen(ex):
    print(line)

