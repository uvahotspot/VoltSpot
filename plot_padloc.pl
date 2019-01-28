#!/usr/bin/perl
use strict;

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
    if ($cur_line =~ /^[VG]\s+?(\d+)\s+?(\d+)/)
    {
        if($1 > $grid_x){
            $grid_x = $1;
        }
        if($2 > $grid_y){
            $grid_y = $2;
        }
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

# file pointers
open $FH, '<', $ARGV[0] or die "Failed with input file '$ARGV[0]': $!";
open my $F_vdd, '>', 'vdd.dat' or die "Failed with output file 'vdd.dat': $!";
open my $F_gnd, '>', 'gnd.dat' or die "Failed with input file 'gnd.dat': $!";

while(<$FH>)
{
    $cur_line = $_;
    if ($cur_line =~ /^V\s+?(\d+)\s+?(\d+)/)
    {
        print $F_vdd "$1\t$2\n";
    }
    elsif ($cur_line =~ /^G\s+?(\d+)\s+?(\d+)/)
    {
        print $F_gnd "$1\t$2\n";
    }
}

$ARGV[0] =~ /^(.+?)\.(.+?)/;
my $VDD_OUT = "$1.vddloc.pdf";
my $GND_OUT = "$1.gndloc.pdf";

open my $F_plot_v, '>', "padloc_v.gpi" or die "Failed with output file 'padloc_v.gpi': $!";

print $F_plot_v "set term postscript enhanced \"Arial\" 20\n";
print $F_plot_v "set pointsize 0.5\n";
print $F_plot_v "set title \"VDD pads location\" font \"Arial,40\"\n";
print $F_plot_v "set size $x_ratio, $y_ratio\n";
print $F_plot_v "set xrange [0:$grid_x]\n";
print $F_plot_v "set yrange [0:$grid_y]\n";
print $F_plot_v "set output \"temp.ps\"\n";
print $F_plot_v "plot \"./vdd.dat\" using 1:2 notitle with points\n";
close $F_plot_v;

system("gnuplot padloc_v.gpi");
system ("ps2pdf temp.ps $VDD_OUT");
system ("rm temp.ps -f");
system ("rm vdd.dat -f");
system ("rm padloc_v.gpi -f");

open my $F_plot_g, '>', "padloc_g.gpi" or die "Failed with output file 'padloc_g.gpi': $!";

print $F_plot_g "set term postscript enhanced \"Arial\" 20\n";
print $F_plot_g "set pointsize 0.5\n";
print $F_plot_g "set title \"GND pads location\" font \"Arial,40\"\n";
print $F_plot_g "set size $x_ratio, $y_ratio\n";
print $F_plot_g "set xrange [0:$grid_x]\n";
print $F_plot_g "set yrange [0:$grid_y]\n";
print $F_plot_g "set output \"temp.ps\"\n";
print $F_plot_g "plot \"./gnd.dat\" using 1:2 notitle with points\n";
close $F_plot_g;

system("gnuplot padloc_g.gpi");
system ("ps2pdf temp.ps $GND_OUT");
system ("rm temp.ps -f");
system ("rm gnd.dat -f");
system ("rm padloc_g.gpi -f");

close $FH;
close $F_vdd;
close $F_gnd;
close $F_plot_v;
close $F_plot_g;
