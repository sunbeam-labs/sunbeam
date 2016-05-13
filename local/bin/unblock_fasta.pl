#!/usr/bin/perl
#
# Author: Christel Chehoud
#
# file: unblock_fasta.pl
# unblocks fasta format files into this format:
# > <ID> 
# <sequencesequencesequencesequencesequencesequencesequencesequencesequence>
# instead of 
# > <ID> 
# <sequence> 
# <sequence>
# <sequence>


$/ = '>';  # read FASTA-delimited records

while (<>) {
    chomp;
    my ($name,@sequence) = split "\n";  # split record into the >id line and several seq lines
  next unless $name =~ /^(\S+)/;      # look for the id at the beginning
  print ">",$name,"\n",@sequence,"\n";       # print  > id, newline, sequence lines and newline
}
