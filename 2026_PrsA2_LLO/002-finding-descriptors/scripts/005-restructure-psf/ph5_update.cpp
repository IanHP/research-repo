# Load your topology
parm /ihome/jdurrant/amm503/ix/durrantlab/Lm-PrsA2-LLO/study/data/004-lm-prsa2-llo-ph5/simulations/01-prep/mol.prmtop

# Change the chain ID of all residues to 'A'
change oresnums :1-273 min 1 max 273
change oresnums :274-546 min 274 max 546
change oresnums :547-1051 min 547 max 1051
change chainid :* A

strip :WAT,Na+,Cl- parmout prsa2-llo-ph5-r1_op.psf

run
