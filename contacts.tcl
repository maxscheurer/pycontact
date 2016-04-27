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
	global sasas
	global sasa_on
	global fnumber
	set fnumber [expr $fnumber + 1]

	set selection1 [atomselect top $sel1]
	set selection2 [atomselect top $sel2]
	set all [atomselect top "$sel1 or $sel2"]
	set really_all [atomselect top "all"]
	set bb_list [$really_all get backbone]
	set resname_list [$really_all get resname]
	set resid_list [$really_all get resid]

	set contactList [measure contacts $cutoff $selection1 $selection2]
	set list1 [lindex $contactList 0]
	set list2 [lindex $contactList 1]

	# set current_results []
	#puts [llength $list1]
	foreach a1 $list1 a2 $list2 {
		set length [measure bond [list $a1 $a2] frame $frame]
		set weight [weight_distance $length]
		# lappend current_results [list $a1 $a2 $length $weight]
		# if {[info exists scorerA($a1)]} {
		# 	set current_a1weight $scorerA($a1)
		# } else {
		# 	set current_a1weight 0.0
		# }

		# if {[info exists scorerB($a2)]} {
		# 	set current_a2weight $scorerB($a2)
		# } else {
		# 	set current_a2weight 0.0
		# }
		# set scorerA($a1) [expr $weight + $current_a1weight]
		# set scorerB($a2) [expr $weight + $current_a2weight]

		#massive speed problems here...
		#set res1 [atomselect top "index $a1"]
		#set res2 [atomselect top "index $a2"]

		#set res1backbone [lindex [$res1 get backbone] 0]
		#set res2backbone [lindex [$res2 get backbone] 0]
		set res1backbone [lindex $bb_list $a1]
		set res2backbone [lindex $bb_list $a2]

		set res1backboneScore 0
		set res1sidechainScore 0

		if {$res1backbone != 0} {
			set res1backboneScore $weight
		} else {
			set res1sidechainScore $weight
		}

		set res2backboneScore 0
		set res2sidechainScore 0

		if {$res2backbone != 0} {
			set res2backboneScore $weight
		} else {
			set res2sidechainScore $weight
		}

		set ss_bb_list [list $res1backboneScore $res1sidechainScore $res2backboneScore $res2sidechainScore]

		#!!! here, one could define more general selections, e.g. for protein-membrane interactions !!!
		# set r1 [$res1 get resid]
		# set n1 [$res1 get resname]
		# set r2 [$res2 get resid]
		# set n2 [$res2 get resname]

		set r1 [lindex $resid_list $a1]
		set n1 [lindex $resname_list $a1]
		set r2 [lindex $resid_list $a2]
		set n2 [lindex $resname_list $a2]
		set key "$n1 $r1 $n2 $r2"

		#SASA
		#more massive speed problems
		if {![info exists sasas($key,$frame)] && $sasa_on} {
			# puts "$key, $frame"
			set residues [atomselect top "$sel1 and resid $r1 or $sel2 and resid $r2"]
			set contact_sasa [measure sasa 1.4 $all -restrict $residues]

			set singleA [atomselect top "$sel1 and resid $r1"]
			set singleB [atomselect top "$sel2 and resid $r2"]

			set singleSasaA [measure sasa 1.4 $selection1 -restrict $singleA]
			set singleSasaB [measure sasa 1.4 $selection2 -restrict $singleB]
			set sasas($key,$frame) [list $contact_sasa [expr $singleSasaA + $singleSasaB]]
			$residues delete
			$singleA delete
			$singleB delete
		}

		# $res1 delete
		# $res2 delete

		if {[info exists residScores($key)]} {
			set current_res_weights $residScores($key)
		} else {
			set current_res_weights []
		}

		lappend current_res_weights [list $frame $weight $ss_bb_list]
		set residScores($key) $current_res_weights
	}
	# lappend results [list $frame $current_results]
	$selection1 delete
	$selection2 delete
	$all delete
	$really_all delete
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
array set sasas {}

set cutoff 5.0
set sasa_on 0
puts "starting contacts"
set fnumber 0
bigdcd contacts auto $traj.dcd
vwait bigdcd_running

# $groupA set beta -5.0
# foreach key [array names scorerA] {
# 	[atomselect top "index $key"] set beta [expr $scorerA($key)/2500.0]
# }
# $groupA writepdb groupA_beta.pdb

# foreach key [array names scorerB] {
# 	[atomselect top "index $key"] set beta [expr $scorerB($key)/2500.0]
# }
# $groupB writepdb groupB_beta.pdb

# set groupA_idx [[atomselect $mol "$sel1"] get index]
# set sumA 0
# foreach idxA groupA_idx {
	
# }

# set outputB [open "outB.dat" w]
# set groupB_res [lsort -unique -integer [[atomselect $mol "$sel2"] get resid]]


# foreach resB $groupB_res {
# 	set ressel [atomselect $mol "resid $resB and $sel2"]
# 	set atoms [$ressel get index]
# 	set sumB 0
# 	foreach idx $atoms {
# 		if {[info exists scorerB($idx)]} {
# 			set current $scorerB($idx)
# 		} else {
# 			set current 0.0
# 		}
# 		set sumB [expr $sumB + $current]
# 	}	
# 	puts $outputB "$resB $sumB"
# }
# close $outputB

set nframes $fnumber
set pyout [open "pyout.dat" w]
puts [array get sasas]
foreach key [array names residScores] {
	set list $residScores($key)
	set keystring ""
	set framestring ""
	append keystring "$key "

	set res1backboneScore 0
	set res1sidechainScore 0
	set res2backboneScore 0
	set res2sidechainScore 0

	for {set i 0} {$i < $nframes} {incr i} {
		set search [lsearch -all -index 0 $list $i]
		if {$search == -1} {
			append framestring "0 "
		} else {
			set framesum 0
			foreach idc $search {
				set framesum [expr $framesum + [lindex $list $idc 1]]
				set res1backboneScore [expr $res1backboneScore + [lindex [lindex $list $idc 2] 0]] ;# maybe change, seems to look complicated
				set res1sidechainScore [expr $res1sidechainScore + [lindex [lindex $list $idc 2] 1]]
				set res2backboneScore [expr $res2backboneScore + [lindex [lindex $list $idc 2] 2]]
				set res2sidechainScore [expr $res2sidechainScore + [lindex [lindex $list $idc 2] 3]]
			}
			append framestring "$framesum "
		}
	}
	append keystring "$res1backboneScore $res1sidechainScore $res2backboneScore $res2sidechainScore "
	append keystring $framestring
	puts $pyout $keystring
}
close $pyout
#set output [open "contacts_$cutoff.dat" w]
#puts $output $results
#close $output

quit
