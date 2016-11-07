set global1 "(protein and segname RN11)"
set global2 "protein"
set restraints "(segname RN11 and same residue as within 5 of segname UBQ)"
set outfile "rn11.dat" 

proc calcSasaPerResidue {globaltext restricttext} {
	#set globaltext "protein"
	#set restricttext "resid 1 2 3"
	set global [atomselect top $globaltext]
	set restrict [atomselect top $restricttext]
	set textprev ""

	array set residues ""
	set resLis ""
	set frames [molinfo top get numframes]
	for {set i 0} {$i < $frames} {incr i} {
                puts "Calculating frame $i of $frames"
		$global frame $i
		$restrict frame $i
		foreach resid [$restrict get resid] resname [$restrict get resname] chain [$restrict get chain] {
			set text "$resname$resid$chain"
			if {$text != $textprev} {
				if {$i == 0} {
					lappend resLis $text
				}
				set selaux [atomselect top "resid $resid and resname $resname and chain $chain" frame $i]
				
				#set aux $residues($text)
				set residues($text,$i) [measure sasa 1.4 $global -restrict $selaux -samples 50]
			
				$selaux delete
			}
			set textprev $text
		}

		set residues(total,$i) [measure sasa 1.4 $global -restrict $restrict -samples 50]

	}

	$global delete
	$restrict delete

	return "{[array get residues]} {$resLis}"


}



set file [open $outfile w+ ]
set aux [calcSasaPerResidue $global1 $restraints]
array set residues1 [lindex $aux 0]
puts [lindex $aux 1]
set resList [lindex $aux 1]

array set residues2 [lindex [calcSasaPerResidue  $global2 $restraints] 0]
set frames [molinfo top get numframes]

#array set resdif "" 
for {set nf 0} {$nf < $frames} {incr nf} {
	for {set i 0} {$i < [llength $resList]} {incr i} {
		puts "$residues1([lindex $resList $i],$nf) - $residues2([lindex $resList $i],$nf)"
		puts $file "[lindex $resList $i],$nf [expr $residues1([lindex $resList $i],$nf) - $residues2([lindex $resList $i],$nf)]"
	}
	puts $file "total,$nf [expr $residues1(total,$nf) - $residues2(total,$nf)]"
puts $file "------------------------------------------------------------------------------"
}
close $file

