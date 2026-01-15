
mol new "D:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run1\\init_coords\\tog_step3_input.psf"
mol addfile "D:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run1\\init_coords\\tog_step5.dcd" step 1 waitfor all
mol addfile "D:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run2\\init_coords\\tog_step5.dcd" step 1 waitfor all
mol addfile "D:\\Research Lab\\MD Sim\\completed_md\\WT_v11\\Run3\\init_coords\\tog_step5.dcd" step 1 waitfor all
#hb
set sel1 [atomselect top "resname SAH"]
set sel2 [atomselect top "protein"]
package require hbonds
hbonds -sel1 $sel1 -sel2 $sel2  -writefile yes -upsel yes -frames all -dist 3.5 -ang 30 -plot no -outdir "D:\\Research Lab\\MD Sim\\completed_md\\data_storage\\WT_G196F\\SAH_h_bonds\\WT" -log allhb.dat -writefile yes -outfile hbonds.dat -polar yes -DA both -type unique -detailout allhb_detailed35_30.dat
exit