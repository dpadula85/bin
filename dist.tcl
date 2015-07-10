# A script for VMD allowing to measure the distance between two residues
# in a protein

proc measure_dist_cent_mass { selstring1 selstring2 {molID {top}} {weight
{}}} {

set sel1 [atomselect $molID $selstring1]
set sel2 [atomselect $molID $selstring2]

if {![llength $weight]} {
set cent_sel1 [measure center $sel1]
set cent_sel2 [measure center $sel2]
} else {
set cent_sel1 [measure center $sel1 weight $weight]
set cent_sel2 [measure center $sel2 weight $weight]
}

set dist [veclength [vecsub $cent_sel1 $cent_sel2]]

$sel1 delete
$sel2 delete
return $dist

}
