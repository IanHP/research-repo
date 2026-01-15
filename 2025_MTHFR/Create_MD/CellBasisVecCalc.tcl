mol new "F:\\Research Lab\\MD Sim\\MD_setup\\WT_v10.1\\initial\\ionized.psf"
mol addfile "F:\\Research Lab\\MD Sim\\MD_setup\\WT_v10.1\\initial\\ionized.pdb"

proc get_box {{molid top}} {

	set all [atomselect $molid all]
	set outfile [open "F:\\Research Lab\\MD Sim\\MD_setup\\WT_v10.1\\initial\\box.txt" w]
	set minmax [measure minmax $all]
	set vec [vecsub [lindex $minmax 1] [lindex $minmax 0]]

	puts "cellBasisVector1 [format "%.4f" [lindex $vec 0]] 0 0"
	puts "cellBasisVector2 0 [format "%.4f" [lindex $vec 1]] 0"
	puts "cellBasisVector3 0 0 [format "%.4f" [lindex $vec 2]]"

	set center [measure center $all]

	puts "cellOrigin [format "%.4f %.4f %.4f" [lindex $center 0] [lindex $center 1] [lindex $center 2]]"

	puts $outfile "cellBasisVector1 [format "%.4f" [lindex $vec 0]] 0 0"
	puts $outfile "cellBasisVector2 0 [format "%.4f" [lindex $vec 1]] 0"
	puts $outfile "cellBasisVector3 0 0 [format "%.4f" [lindex $vec 2]]"
	puts $outfile "cellOrigin [format "%.4f %.4f %.4f" [lindex $center 0] [lindex $center 1] [lindex $center 2]]"

	$all delete
}

get_box

exit

