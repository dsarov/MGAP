#! /usr/bin/perl

my $length;
my $totalLength;
my @arr;
while(){
   chomp;
   if(/>/){
   push (@arr, $length);
   $totalLength += $length;
   $length=0;
   next;
  }
  $length += length($_);
}

close(FH);

my @sort = sort {$b <=> $a} @arr;
my $n50;
my $L50;
foreach my $val(@sort){
     $n50+=$val;
     $L50++;
      if($n50 >= $totalLength/2){
         print "N50 length is $n50 and N50 value is: $val and L50 is $L50\n";
         last;
     }
}
