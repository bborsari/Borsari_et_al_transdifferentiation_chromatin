#!/usr/bin/perl -w

# Author: Ramil Nurtdinov

use strict;

my $need_logarithm       = 0;    #Is it nessesary to make log2(1+FPKM)? 0 - no, 1 - yes
my $fold_change          = 1;    #Based on log2 scale
my $peak_heght           = 1;    #Heght of the peak in log2(1+FPKM) scale

my $prof_threshold       = 0.50; #Maximum proportion of gene interval that can be treated as a deviation
my $comp_threshold       = 0.25; #Maximum amount of non compensated deviation
my $low_expression_mean  = 0.5;  #The value of mean expression to consider gene low (none) expressed, FPRM = 0.414213562
my $high_expression_mean = 4;    #The value of mean expression to consider gene high expressed, FPRM = 15

my $time_series_lenght = 12;   #How many timepoints we have in series
my $time_series_last   = $time_series_lenght - 1; #Array starts from 0

my $file_name = $ARGV[0];
unless($file_name)
{
    print "Please provide a name of file with FPKM or log(FPKM).\nIn case of log(FPKM) please note that algorithm may\nnot work good with negative values to avoid\nthis I am using pseudocounts log2(1+FPKM)\n";
    exit;
}

open FILE_IN,"$file_name" or die "cannot open file $file_name\n";
my $string = <FILE_IN>;
chomp $string;
$string=~s/\r//;
my @Header = split "\t", $string;
shift @Header;

open FILE_OUT,">$file_name.classification";
print FILE_OUT "gene_id\tclass\ttimepoint\n";

my %Classes = ();
while($string = <FILE_IN>)
{
    chomp $string;$string=~s/\r//;
    my @Data = split "\t", $string;
    my $gene_id = shift @Data;
    if($need_logarithm==1)
    {
	@Data = map { &get_log($_) } @Data[0..$time_series_last];
    }
    @Data = @Data[0..$time_series_last];
    my ($class,$extra) = &classify_timeseries(\@Data);
    $Classes{$class}++;

    print FILE_OUT "$gene_id\t$class\t$extra\n";
}

map { print $_,"\t",$Classes{$_},"\n" } ("low expression","moderate expression","high expression");
print "\n";

map { print $_,"\t",$Classes{$_},"\n" } ("upregulation","downregulation","peaking","bending");
print "\n";

map { print $_,"\t",$Classes{$_},"\n" } ("variable");


sub classify_timeseries
{
    my $array_ref = $_[0];
    my @TimeSeries = map { "$_" } @{$array_ref};
    my @Sorted     = sort { $b <=> $a } map { "$_" } @TimeSeries;
    my $interval   = $Sorted[0] - $Sorted[-1];

    my $sum = 0;
    map { $sum+=$_ } @Sorted;
    $sum/=$time_series_lenght;
    return ("low expression", "flat") if $sum < $low_expression_mean;

    if($interval < $fold_change)
    {
	return ("moderate expression","flat") if $sum < $high_expression_mean;
	return ("high expression","flat");
    }

    my ($upregulation,$downregulation) = time_series($array_ref);

    return ("upregulation","NA")   if $upregulation > $downregulation and $downregulation <= $prof_threshold and $upregulation - $downregulation >= (1-$comp_threshold);
    return ("downregulation","NA") if $downregulation > $upregulation and $upregulation   <= $prof_threshold and $downregulation - $upregulation >= (1-$comp_threshold);

    my ($type,$timepoint) = &get_peaking($array_ref);
    return ($type,$timepoint) if $type eq "peaking";
    
    ($type,$timepoint) = &get_bending($array_ref);
    return ($type,$timepoint) if $type eq "bending";

    return ("variable","variable");
}

sub get_peaking
{
    my $array_ref  = $_[0];
    my @TimeSeries = map { "$_" } @{$array_ref};

    my @Sorted = sort { $b->[0] <=> $a->[0] } map { ["$TimeSeries[$_]","$_"] } (0..$time_series_last);
    my ($peaking_value,$peaking_point) = @{$Sorted[0]};

    return ("bad","NA") if $peaking_point < 1 or $peaking_point > $time_series_last-1;

    my $delta_begin = $peaking_value - $TimeSeries[0];
    my $delta_end   = $peaking_value - $TimeSeries[$time_series_last];

    return ("bad","NA") if $delta_begin < $peak_heght or $delta_end < $peak_heght;


    my @New_values=();
    foreach my $index (0..$peaking_point)
    {
        push @New_values,$TimeSeries[$index];
    }
    $peaking_point++;
    foreach my $index ($peaking_point..$time_series_last)
    {
        my $new_exp = $peaking_value + $peaking_value - $TimeSeries[$index];
        push @New_values,sprintf("%.5f",$new_exp);
    }

    my ($upregulation,$downregulation) = &time_series(\@New_values);
    return ("peaking",$Header[$peaking_point-1]) if $upregulation > $downregulation and $downregulation <= $prof_threshold and $upregulation - $downregulation >= (1-$comp_threshold);
    return ("bad","NA");
}


sub get_bending
{
    my $array_ref  = $_[0];
    my @TimeSeries = map { "$_" } @{$array_ref};

    my @Sorted = sort { $a->[0] <=> $b->[0] } map { ["$TimeSeries[$_]","$_"] } (0..$time_series_last);
    my ($bending_value,$bending_point) = @{$Sorted[0]};

    return "bad" if $bending_point < 2 or $bending_point > $time_series_last-1;

    my $delta_begin = $TimeSeries[0] - $bending_value;
    my $delta_end   = $TimeSeries[9] - $bending_value;

    return "bad" if $delta_begin < $peak_heght or $delta_end < $peak_heght;


    my @New_values=();
    foreach my $index (0..$bending_point)
    {
        push @New_values,$TimeSeries[$index];
    }
    $bending_point++;
    foreach my $index ($bending_point..$time_series_last)
    {
        my $new_exp = $bending_value + $bending_value - $TimeSeries[$index];
        push @New_values,sprintf("%.5f",$new_exp);
    }

    my ($upregulation,$downregulation) = &time_series(\@New_values);
    return ("bending",$Header[$bending_point-1]) if $downregulation > $upregulation and $upregulation   <= $prof_threshold and $downregulation - $upregulation >= (1-$comp_threshold);
    return ("bad","NA");
}



sub time_series
{
    my $array_ref  = $_[0];
    my @TimeSeries = map { "$_" } @{$array_ref};
    my @Sorted     = sort { $b <=> $a } map { "$_" } @TimeSeries;

    my $interval   = $Sorted[0] - $Sorted[-1];
    my ($upregulation,$downregulation) = (0,0);

    foreach my $index (1..$time_series_last)
    {
        if($TimeSeries[$index] >= $TimeSeries[$index-1])
        {
            $upregulation   += $TimeSeries[$index] - $TimeSeries[$index-1];
        }
        else
        {
            $downregulation += $TimeSeries[$index-1] - $TimeSeries[$index];
        }
    }
    $upregulation   = sprintf("%.4f",$upregulation/$interval);
    $downregulation = sprintf("%.4f",$downregulation/$interval);
    return ($upregulation,$downregulation);
}




sub get_log
{
    my $expression = $_[0];
    my $new_expression = log(1+$expression)/log(2);
    return sprintf("%.6f",$new_expression) if $new_expression < 10;
    return sprintf("%.5f",$new_expression);
}


