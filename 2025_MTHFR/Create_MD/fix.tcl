mol load psf "D:\\Research Lab\\MD Sim\\MD_setup\\WT_v7\\initial\\ionized.psf" pdb "D:\\Research Lab\\MD Sim\\MD_setup\\WT_v7\\initial\\ionized.pdb"

set s1 [atomselect top all]
$s1 set beta 0
$s1 set occupancy 0

set s2 [atomselect top protein]
$s2 set beta 1

$s1 writepdb "D:\\Research Lab\\MD Sim\\MD_setup\\WT_v7\\reffix.pdb"