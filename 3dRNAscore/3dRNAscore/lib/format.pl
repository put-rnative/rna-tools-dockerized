#! /usr/bin/perl

if ($ARGV[0] eq "-check") {
	$check = 1;
	open FILE, "< $ARGV[1]";
} else {
	$check = 0;
	open FILE, "< $ARGV[0]";
}

$res_number = 1;
while (<FILE>) {
	last if /ENDMDL/;
	chomp;
	next if not /^ATOM/;
	$line = $_;
	$line =~ s/'/*/g;

	$tmp = substr($line, 16, 1);
	next if not $tmp eq 'A' and not $tmp eq ' ';

	$residue = &del_space(uc(substr($line, 17, 3)));
	if ($residue eq "A" or $residue eq "RA" or $residue eq "RA5" or $residue eq "RA3" or $residue eq "A5" or $residue eq "A3") {
		$residue = "A";
	} elsif ($residue eq "U" or $residue eq "RU" or $residue eq "RU5" or $residue eq "RU3" or $residue eq "U5" or $residue eq "U3") {
		$residue = "U";
	} elsif ($residue eq "G" or $residue eq "RG" or $residue eq "RG5" or $residue eq "RG3" or $residue eq "G5" or $residue eq "G3") {
		$residue = "G"
	} elsif ($residue eq "C" or $residue eq "RC" or $residue eq "RC5" or $residue eq "RC3" or $residue eq "C5" or $residue eq "C3") {
		$residue = "C"
	} else {
		next;
	}
	
	$n++;

	$chain = &del_space(substr($line, 21, 1));
	if ($chain eq "") {
		$chain = "X";
	}

	$residue_num = &del_space(substr($line, 22, 5));
	if ((not $residue_num eq $old_residue_num or not $chain eq $old_chain or not $residue eq $old_residue) and $n != 1) {
		if (&print_residue) {
			undef %x;
			undef %y;
			undef %z;
			undef @print_list;			
		};
		if (not $chain eq $old_chain) {
			$echo_p = 0;
		} else {
			$echo_p = 1;
		}
	}
	$old_residue = $residue;
	$old_residue_num = $residue_num;

	$atom = &del_space(substr($line, 12, 4));
	$atom = "O1P" if $atom eq "OP1";
	$atom = "O2P" if $atom eq "OP2";
	if (defined($x{$atom})) {
		die "too many atoms: $atom!\n";
	}
	$x{$atom} = substr($line, 30, 8);
	$y{$atom} = substr($line, 38, 8);
	$z{$atom} = substr($line, 46, 8);
	
	print "TER\n" if not $chain eq $old_chain and $n != 1 and not $check;
	$old_chain = $chain;
}
&print_residue;
print "TER\n" if not $check;
close FILE;

sub del_space {
	my $a = $_[0];
	$a =~ s/^\s+//g;
	$a =~ s/\s+$//g;
	$a;
}

sub print_residue {
	$res_n++;
	# check residue
	for $j ("O5*", "C5*", "C4*", "O4*", "C3*", "O3*", "C2*", "O2*", "C1*") {
		if (not defined($x{$j})) {
			print "chain ".$old_chain." residue".$old_residue_num.' '.$old_residue.' '."no atom: $j\n" if $check;
			$res_n--;
			return 1;
		}
		$print_list[@print_list] = $j;
	}
	if ($old_residue eq "A") {
		for $j ("N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4") {
			if (not defined($x{$j})) {
				print "chain ".$old_chain." residue".$old_residue_num.' '.$old_residue.' '."no atom: $j\n" if $check;
				$res_n--;
				return 1;
			}
			$print_list[@print_list] = $j;
		}
	} elsif ($old_residue eq "U") {
		for $j ("N1", "C2", "N3", "O2", "C4", "O4", "C5", "C6") {
			if (not defined($x{$j})) {
				print "chain ".$old_chain." residue".$old_residue_num.' '.$old_residue.' '."no atom: $j\n" if $check;
				$res_n--;
				return 1;
			}
			$print_list[@print_list] = $j;
		}
	} elsif ($old_residue eq "G") {
		for $j ("N9", "C8", "N7", "C5", "C6", "O6", "N1", "C2", "N2", "N3", "C4") {
			if (not defined($x{$j})) {
				print "chain ".$old_chain." residue".$old_residue_num.' '.$old_residue.' '."no atom: $j\n" if $check;
				$res_n--;
				return 1;
			}
			$print_list[@print_list] = $j;
		}
	} elsif ($old_residue eq "C") {
		for $j ("N1", "C2", "N3", "O2", "C4", "N4", "C5", "C6") {
			if (not defined($x{$j})) {
				print "chain ".$old_chain." residue".$old_residue_num.' '.$old_residue.' '."no atom: $j\n" if $check;
				$res_n--;
				return 1;
			}
			$print_list[@print_list] = $j;
		}
	}

	# print residue
	if (not $check) {
		if (defined($x{"P"}) and defined($x{"O1P"}) and defined($x{"O2P"}) and $res_n != 1 and $echo_p == 1) {
			printf "ATOM %6d  %-3s%4s %1s%4d    %8.3f%8.3f%8.3f\n", ++$atom_num, "P", $old_residue, $old_chain, $res_number, $x{"P"}, $y{"P"}, $z{"P"};
			printf "ATOM %6d  %-3s%4s %1s%4d    %8.3f%8.3f%8.3f\n", ++$atom_num, "O1P", $old_residue, $old_chain, $res_number, $x{"O1P"}, $y{"O1P"}, $z{"O1P"};
			printf "ATOM %6d  %-3s%4s %1s%4d    %8.3f%8.3f%8.3f\n", ++$atom_num, "O2P", $old_residue, $old_chain, $res_number, $x{"O2P"}, $y{"O2P"}, $z{"O2P"};
		}
		for $i (@print_list) {
			printf "ATOM %6d  %-3s%4s %1s%4d    %8.3f%8.3f%8.3f\n", ++$atom_num, $i, $old_residue, $old_chain, $res_number, $x{$i}, $y{$i}, $z{$i};
		}
		$res_number++;
	}

	undef %x;
	undef %y;
	undef %z;
	undef @print_list;
}





