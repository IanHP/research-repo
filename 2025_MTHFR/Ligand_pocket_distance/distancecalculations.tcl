# Contact for residue
# Sema Zeynep Yılmaz 26.10.22
#PR65 başlangıç yapısı
#MolID 0
# CA distance measure script
# between two residue through trajectory

proc int_dist {sele1 sele2 outfile} {
            set sel1a [atomselect top "$sele1"]
            set sel2a [atomselect top "$sele2"]
			
set log [open $outfile w]

    set num [molinfo top get numframes]

    for {set i 0} {$i < $num} {incr i} {
        puts "frame $i of $num"
        $sel1a frame $i
        $sel2a frame $i
		set sel1 [measure center $sel1a weight mass]
		set sel2 [measure center $sel2a weight mass]

        set distance [veclength [vecsub $sel1 $sel2]]
        puts $log "$i\t$distance"
    }
    close $log
	}


# WT
mol new "F:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run1\\init_coords\\tog_step3_input.psf"
mol addfile "F:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run1\\init_coords\\tog_step5.dcd" waitfor all
int_dist "protein and chain A and resid 428 to 440 462 to 492 506 to 517 573 to 598 604 to 623 and name CA " "resname SAH and noh" posedistance_WT1.dat
mol delete all

mol new "F:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run2\\init_coords\\tog_step3_input.psf"
mol addfile "F:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run2\\init_coords\\tog_step5.dcd" waitfor all
int_dist "protein and chain A and resid 428 to 440 462 to 492 506 to 517 573 to 598 604 to 623 and name CA " "resname SAH and noh" posedistance_WT2.dat
mol delete all

mol new "F:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run3\\init_coords\\tog_step3_input.psf"
mol addfile "F:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run3\\init_coords\\tog_step5.dcd" waitfor all
int_dist "protein and chain A and resid 428 to 440 462 to 492 506 to 517 573 to 598 604 to 623 and name CA " "resname SAH and noh" posedistance_WT3.dat
mol delete all


# G196F
mol new "F:\\Research Lab\\MD Sim\\completed_md\\G196F_v1\\Run1\\init_coords\\tog_step3_input.psf"
mol addfile "F:\\Research Lab\\MD Sim\\completed_md\\G196F_v1\\Run1\\init_coords\\tog_step5.dcd" waitfor all
int_dist "protein and chain A and resid 428 to 440 462 to 492 506 to 517 573 to 598 604 to 623 and name CA " "resname SAH and noh" posedistance_G196F1.dat
mol delete all

mol new "F:\\Research Lab\\MD Sim\\completed_md\\G196F_v1\\Run2\\init_coords\\tog_step3_input.psf"
mol addfile "F:\\Research Lab\\MD Sim\\completed_md\\G196F_v1\\Run2\\init_coords\\tog_step5.dcd" waitfor all
int_dist "protein and chain A and resid 428 to 440 462 to 492 506 to 517 573 to 598 604 to 623 and name CA " "resname SAH and noh" posedistance_G196F2.dat
mol delete all

mol new "F:\\Research Lab\\MD Sim\\completed_md\\G196F_v1\\Run3\\init_coords\\tog_step3_input.psf"
mol addfile "F:\\Research Lab\\MD Sim\\completed_md\\G196F_v1\\Run3\\init_coords\\tog_step5.dcd" waitfor all
int_dist "protein and chain A and resid 428 to 440 462 to 492 506 to 517 573 to 598 604 to 623 and name CA " "resname SAH and noh" posedistance_G196F3.dat
mol delete all



exit

 
