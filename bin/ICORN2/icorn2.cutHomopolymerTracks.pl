use strict;

### script will cut all the reads that have at last 13 T or A
while (<STDIN>){
	my $r=$_;
	$_=<STDIN>;
	if (/A{13,}/ || /T{13,}/){
			my $pos=@-[0];
		if ($pos > 50){
			print $r;
			print substr($_,0,($pos+5))."\n";
			$_=<STDIN>; print ;	
			$_=<STDIN>;print substr($_,0,($pos+5))."\n";
		}
		else {
				print $r;
				print "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n+\n"	;
			print "#########################################################\n";
			
			$_=<STDIN>;$_=<STDIN>;
		}
	}
	else {
		print $r;
		print ;
		$_=<STDIN>; print ;	
		$_=<STDIN>; print ;	
	}
}
