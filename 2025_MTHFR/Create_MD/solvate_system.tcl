set midlist {}

set mol [mol new "F:\\Research Lab\\MD Sim\\MD_setup\\WT_v10.1\\initial\\autopsf.psf" waitfor all]
mol addfile "F:\\Research Lab\\MD Sim\\MD_setup\\WT_v10.1\\initial\\autopsf.pdb"
lappend midlist $mol

set mol [mol new "F:\\Research Lab\\MD Sim\\MD_setup\\WT_v10.1\\ligands_toppar\\fad.psf" waitfor all]
mol addfile "F:\\Research Lab\\MD Sim\\MD_setup\\WT_v10.1\\ligands_toppar\\fad.pdb" $mol
lappend midlist $mol

set mol [mol new "F:\\Research Lab\\MD Sim\\MD_setup\\WT_v10.1\\ligands_toppar\\sah.psf" waitfor all]
mol addfile "F:\\Research Lab\\MD Sim\\MD_setup\\WT_v10.1\\ligands_toppar\\sah.pdb" $mol
lappend midlist $mol

# do the magic
set mol [::TopoTools::mergemols $midlist]

set a [atomselect top all]

set c [measure center $a weight mass]

$a moveby [vecsub {0 0 0} $c]

set a [atomselect top all]

$a writepsf aligned.psf
$a writepdb aligned.pdb

package require solvate

solvate aligned.psf aligned.pdb -o solvate -s WT -minmax {{-69.5 -69.5 -69.5} {69.5 69.5 69.5}} -x 0 -y 20 -z -10 +x 0 +y -20 +z 10 -b 2.4