
### WT - RUN 1
mol new "F:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run1\\init_coords\\tog_step3_input.psf"
mol addfile "F:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run1\\init_coords\\tog_step3_full_FAD_align.dcd" waitfor all

set numbers {119 187 237 246 304 318 328 946 1020 1168 1244 1554}
foreach i $numbers {
    set c2 [atomselect top "all" frame $i]
    $c2 writepdb "$i.pdb"
}

mol delete all

mol new 119.pdb
set numbers {187 237 246 304 318 328 946 1020 1168 1244 1554}
foreach i $numbers {
    mol addfile $i.pdb
}

set c2 [atomselect top "all"]
animate write dcd wt1_selected_frame.dcd sel $c2 skip 1 waitfor all




### WT - RUN 2
mol new "F:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run2\\init_coords\\tog_step3_input.psf"
mol addfile "F:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run2\\init_coords\\tog_step3_full_FAD_align.dcd" waitfor all

set numbers {119 187 237 246 304 318 328 946 1020 1168 1244 1554}
foreach i $numbers {
    set c2 [atomselect top "all" frame $i]
    $c2 writepdb "$i.pdb"
}

mol delete all

mol new 119.pdb
set numbers {187 237 246 304 318 328 946 1020 1168 1244 1554}
foreach i $numbers {
    mol addfile $i.pdb
}

set c2 [atomselect top "all"]
animate write dcd wt2_selected_frame.dcd sel $c2 skip 1 waitfor all





### WT - RUN 3
mol new "F:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run3\\init_coords\\tog_step3_input.psf"
mol addfile "F:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run3\\init_coords\\tog_step3_full_FAD_align.dcd" waitfor all

set numbers {119 187 237 246 304 318 328 946 1020 1168 1244 1554}
foreach i $numbers {
    set c2 [atomselect top "all" frame $i]
    $c2 writepdb "$i.pdb"
}

mol delete all

mol new 119.pdb
set numbers {187 237 246 304 318 328 946 1020 1168 1244 1554}
foreach i $numbers {
    mol addfile $i.pdb
}

set c2 [atomselect top "all"]
animate write dcd wt3_selected_frame.dcd sel $c2 skip 1 waitfor all







### MUT - RUN 1
mol new "F:\\Research Lab\\MD Sim\\completed_md\\G196F_v1\\Run1\\init_coords\\tog_step3_input.psf"
mol addfile "F:\\Research Lab\\MD Sim\\completed_md\\G196F_v1\\Run1\\init_coords\\tog_step3_full_FAD_align.dcd" waitfor all

set numbers {119 187 237 246 304 318 328 946 1020 1168 1244 1554}
foreach i $numbers {
    set c2 [atomselect top "all" frame $i]
    $c2 writepdb "$i.pdb"
}

mol delete all

mol new 119.pdb
set numbers {187 237 246 304 318 328 946 1020 1168 1244 1554}
foreach i $numbers {
    mol addfile $i.pdb
}

set c2 [atomselect top "all"]
animate write dcd mut1_selected_frame.dcd sel $c2 skip 1 waitfor all




### MUT - RUN 2
mol new "F:\\Research Lab\\MD Sim\\completed_md\\G196F_v1\\Run2\\init_coords\\tog_step3_input.psf"
mol addfile "F:\\Research Lab\\MD Sim\\completed_md\\G196F_v1\\Run2\\init_coords\\tog_step3_full_FAD_align.dcd" waitfor all

set numbers {119 187 237 246 304 318 328 946 1020 1168 1244 1554}
foreach i $numbers {
    set c2 [atomselect top "all" frame $i]
    $c2 writepdb "$i.pdb"
}

mol delete all

mol new 119.pdb
set numbers {187 237 246 304 318 328 946 1020 1168 1244 1554}
foreach i $numbers {
    mol addfile $i.pdb
}

set c2 [atomselect top "all"]
animate write dcd mut2_selected_frame.dcd sel $c2 skip 1 waitfor all




### MUT - RUN 3
mol new "F:\\Research Lab\\MD Sim\\completed_md\\G196F_v1\\Run3\\init_coords\\tog_step3_input.psf"
mol addfile "F:\\Research Lab\\MD Sim\\completed_md\\G196F_v1\\Run3\\init_coords\\tog_step3_full_FAD_align.dcd" waitfor all

set numbers {119 187 237 246 304 318 328 946 1020 1168 1244 1554}
foreach i $numbers {
    set c2 [atomselect top "all" frame $i]
    $c2 writepdb "$i.pdb"
}

mol delete all

mol new 119.pdb
set numbers {187 237 246 304 318 328 946 1020 1168 1244 1554}
foreach i $numbers {
    mol addfile $i.pdb
}

set c2 [atomselect top "all"]
animate write dcd mut3_selected_frame.dcd sel $c2 skip 1 waitfor all

exit