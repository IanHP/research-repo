mol load psf "F:\\Research Lab\\MD Sim\\MD_setup\\WT_v10.1\\initial\\ionized.psf" pdb "F:\\Research Lab\\MD Sim\\MD_setup\\WT_v10.1\\initial\\ionized.pdb"

set all [atomselect top all]
$all set beta 0
$all set occupancy 0

set A [atomselect top "chain A and backbone and not hydrogen"]
$A set beta 1

set B [atomselect top "chain A and sidechain and not hydrogen"]
$B set beta 0.5

$all writepdb consref.pdb
exit