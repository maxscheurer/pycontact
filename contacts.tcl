# contact analysis script for proteins/proteins and proteins/membranes
# Author: Maximilian Scheurer, mail: mscheurer@ks.uiuc.edu
# April 2016

proc weight_distance {dist} \
{
	return [expr (1.0)/(1.0 + exp(5.0*($dist-4.0)))]
}

proc contacts {frame} \
{
	puts $frame
	global cutoff
	global sel1
	global sel2
	# global results
	global scorerA
	global scorerB
	global residScores

	set selection1 [atomselect top $sel1]
	set selection2 [atomselect top $sel2]

	set contactList [measure contacts $cutoff $selection1 $selection2]
	set list1 [lindex $contactList 0]
	set list2 [lindex $contactList 1]

	# set current_results []
	puts [llength $list1]
	foreach a1 $list1 a2 $list2 {
		set length [measure bond [list $a1 $a2] frame $frame]
		set weight [weight_distance $length]
		# lappend current_results [list $a1 $a2 $length $weight]
		if {[info exists scorerA($a1)]} {
			set current_a1weight $scorerA($a1)
		} else {
			set current_a1weight 0.0
		}

		if {[info exists scorerB($a2)]} {
			set current_a2weight $scorerB($a2)
		} else {
			set current_a2weight 0.0
		}
		set scorerA($a1) [expr $weight + $current_a1weight]
		set scorerB($a2) [expr $weight + $current_a2weight]

		#massive speed problems here...
		set res1 [atomselect top "index $a1"]
		set res2 [atomselect top "index $a2"]

		#!!! here, one could define more general selections, e.g. for protein-membrane interactions !!!
		set r1 [$res1 get resid]
		set n1 [$res1 get resname]
		set r2 [$res2 get resid]
		set n2 [$res2 get resname]
		set key "$n1 $r1 $n2 $r2"

		$res1 delete
		$res2 delete

		if {[info exists residScores($key)]} {
			set current_res_weights $residScores($key)
		} else {
			set current_res_weights []
		}
		lappend current_res_weights [list $frame $weight]
		set residScores($key) $current_res_weights
	}
	# lappend results [list $frame $current_results]
	$selection1 delete
	$selection2 delete
}

source bigdcd.tcl

set psf rpn11_ubq_interface-ionized
set traj short
set mol [mol new $psf.psf type psf waitfor all]

set sel1 "segname RN11 and noh"
set sel2 "segname UBQ and noh"

mol addfile $psf.pdb molid $mol
set groupA [atomselect top "$sel1"]
set groupB [atomselect top "$sel2"]

array set scorerA {}
array set scorerB {}

array set residScores {}

set cutoff 3.0

bigdcd contacts auto $traj.dcd
vwait bigdcd_running

$groupA set beta -5.0
foreach key [array names scorerA] {
	[atomselect top "index $key"] set beta [expr $scorerA($key)/2500.0]
}
$groupA writepdb groupA_beta.pdb

foreach key [array names scorerB] {
	[atomselect top "index $key"] set beta [expr $scorerB($key)/2500.0]
}
$groupB writepdb groupB_beta.pdb

# set groupA_idx [[atomselect $mol "$sel1"] get index]
# set sumA 0
# foreach idxA groupA_idx {
	
# }

set outputB [open "outB.dat" w]
set groupB_res [lsort -unique -integer [[atomselect $mol "$sel2"] get resid]]


foreach resB $groupB_res {
	set ressel [atomselect $mol "resid $resB and $sel2"]
	set atoms [$ressel get index]
	set sumB 0
	foreach idx $atoms {
		if {[info exists scorerB($idx)]} {
			set current $scorerB($idx)
		} else {
			set current 0.0
		}
		set sumB [expr $sumB + $current]
	}	
	puts $outputB "$resB $sumB"
}
close $outputB

set nframes 50
set pyout [open "pyout.dat" w]
foreach key [array names residScores] {
	set list $residScores($key)
	set keystring ""
	append keystring "$key "
	for {set i 0} {$i < $nframes} {incr i} {
		set search [lsearch -all -index 0 $list $i]
		if {$search == -1} {
			append keystring "0 "
		} else {
			set framesum 0
			foreach idc $search {
				set framesum [expr $framesum + [lindex $list $idc 1]]
			}
			append keystring "$framesum "
		}
	}
	puts $pyout $keystring
}
close $pyout
#set output [open "contacts_$cutoff.dat" w]
#puts $output $results
#close $output

quit
