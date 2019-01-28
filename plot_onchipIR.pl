#!/usr/bin/perl
use strict;
use warnings;

my $numArgs;
$numArgs = $#ARGV + 1;
if ($numArgs!=1)
{
    print "wrong number of arguments!\n Usage: this.perl inputfile\n";
    exit
}

# get grid size
open my $FH, '<', $ARGV[0] or die "Failed with input file '$ARGV[0]': $!";
my $cur_line;

my $grid_x = 0;
my $grid_y = 0;

while(<$FH>)
{
    $cur_line = $_;
    if ($cur_line =~ /^(\d+)\s+?(\d+)\s+?([\.\d]+)/)
    {
        if($1 > $grid_x){
            $grid_x = $1;
        }
        if($2 > $grid_y){
            $grid_y = $2;
        }
    }
    elsif (($cur_line =~ /^#/) || ($cur_line =~ /^\s+$/))
    {
        # do nothing
    }
    else
    {
        print "$cur_line\n";
        die "Unrecognized line\n";
    }
}
close $FH;

my $x_ratio = 1;
my $y_ratio = 1;

if($grid_x > $grid_y) {
    $y_ratio = $grid_y/$grid_x;
} else {
    $x_ratio = $grid_x/$grid_y;
}

# generate dat file
open $FH, '<', $ARGV[0] or die "Failed with input file '$ARGV[0]': $!";
open my $F_dat, '>', 'IR.dat' or die "Failed with output file 'IR.dat': $!";

while(<$FH>)
{
    $cur_line = $_;
    if ($cur_line =~ /^(\d+)\s+?(\d+)\s+?([\.\d]+)/)
    {
        print $F_dat "$1\t$2\t$3\n";
    }
}
close $F_dat;
close $FH;

my ($plot_min, $plot_max);
$plot_min = '*'; $plot_max = '*';

$ARGV[0] =~ /^(.+?)\.(.+?)/;
my $out_file = "$1";

open my $F_plot, '>', "temp.gpi" or die "Failed with output file 'temp.gpi': $!";
print $F_plot "set term gif font arial 14\n";
print $F_plot "set style rectangle back fc lt -3 fillstyle  solid 1.00 border -1\n";
print $F_plot "unset key\n";
print $F_plot "set view map\n";
print $F_plot "set xtics border in scale 0,0 mirror norotate  offset character 0, 0, 0\n";
print $F_plot "set ytics border in scale 0,0 mirror norotate  offset character 0, 0, 0\n";
print $F_plot "unset xtics\n";
print $F_plot "unset ytics\n";
print $F_plot "set size $x_ratio, $y_ratio\n";
print $F_plot "set xrange [0 : $grid_x ]\n";
print $F_plot "set yrange [0 : $grid_y ]\n";
print $F_plot "set ylabel  offset character 0, 0, 0 font \"\" textcolor lt -1 rotate by 90\n";
print $F_plot "set cblabel \"IR drop(%Vdd)\"\n";
print $F_plot "set cblabel  offset character 0, 0, 0 font \"\" textcolor lt -1 rotate by 90\n";
print $F_plot "set cbrange [ $plot_min : $plot_max ] noreverse nowriteback\n";
print $F_plot "set palette rgbformulae 21, 22, 23\n";
print $F_plot "set output \"$out_file.gif\"\n";
print $F_plot "plot \"IR.dat\" using 1:2:3 with image\n";
close $F_plot;
system("gnuplot temp.gpi");
system ("rm temp.gpi -f");
close $F_plot;

system ("rm IR.dat -f");
