#!/bin/python
# NOTE: crimm requires seems to require mmcif files for modeling, which need to be downloaded from the PDB
from crimm.Fetchers import ( fetch_rcsb )
#from crimm.IO import PDBParser
from crimm.IO import get_pdb_str ;# to print pdb-formatted lines
from crimm.Modeller.LoopBuilder import ChainLoopBuilder ;
#from Bio.PDB.PDBIO import PDBIO
#from Bio.PDB.PDBIO import Select
from sys import exit

pdb = fetch_rcsb('1hzh')

for chain in pdb.get_chains():
 print("Chain ID:",chain.get_id())
 if chain.chain_type=='Polypeptide(L)':
  print("Full seq:\n>",chain.seq) ;# fails for sugars
  print("Canonical seq:\n>",chain.can_seq)
  print("Missing res:\n>",chain.missing_res)
  print("Masked seq:\n>",chain.masked_seq)
  print(">",end=' '); chain.masked_seq.show();

# create model for chain B

chain = pdb['B'] ;# maybe no scalar key if taking top model
chain.masked_seq.show()
looper = ChainLoopBuilder(chain)
#looper.build_from_homology(max_num_match = 10, identity_score_cutoff = 0.9) ;# fails
looper.build_from_alphafold(include_terminal=False);

full_chain = looper.get_chain()
full_chain.truncate_missing_terminal();
chain=full_chain;
print("Full seq:\n>",chain.seq) ;# fails for sugars
print("Canonical seq:\n>",chain.can_seq)
print("Missing res:\n>",chain.missing_res)
print("Masked seq:\n>",chain.masked_seq)
print(">",end=' '); chain.masked_seq.show();

# write modeled chain coordinate to PDB file
# note that chains, residues and atomd will be renumbered, sigh
# using get_pdb_str:
pdb_str=get_pdb_str(full_chain)
outname="HC2-AF2.pdb";
with open(outname, "w") as of:
  of.write(pdb_str)

exit();# no need to continue

# to write original file with missing loops
# NOTE: it is useless because it renumbers everything, despite having made no changes
pdb_str=get_pdb_str(pdb['B']) 
outname="HC2-MISSING.pdb";
with open(outname, "w") as of:
  of.write(pdb_str)

# using Biopython:
pdbout=PDBIO() ; # init  PDB object
#pdbout.set_structure(pdb) ;# whole thing, renumbers residues
pdbout.set_structure(full_chain) ;# single chain
pdbout.save('HC2-AF2I.pdb',preserve_atom_numbering=True)
# the two writes above are identical, except for the numbering of the TER string at the end
