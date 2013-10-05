#!/usr/bin/perl

use strict;
use warnings;

use POSIX qw(strftime);

use Math::GSL::SF qw/:all/;
use Math::GSL::Randist qw/:all/;

use AnyEvent;
use AnyEvent::ForkManager;

use File::Path;

use Getopt::Std;

my %opts;
getopt('mcepsao',\%opts);

&check_args_user(\%opts);

my ($norm,$rgb,$aN,$scan_out,$track_out);
($norm)                 = 'c';
($rgb,$aN)              = ('255,127,36',"b$norm");
($scan_out,$track_out)  = ("$opts{o}/scan_out_$norm","$opts{o}/tracks_$norm");

my (%events_by_chr,%size_by_chr,%ss_aN,$genome_size,$events_in_genome,@chrs);
&load_chr_size(\$genome_size,\%size_by_chr,$opts{c});
&load_sample(\$events_in_genome,\%events_by_chr,$opts{e});
&calc_ss_by_cg(\%size_by_chr,\%events_by_chr,\%ss_aN,$opts{m},$opts{p},$opts{e},$opts{a},'chr',$scan_out);

@chrs   = (keys %events_by_chr);
&mk_track(\@chrs,\%ss_aN,$opts{e},$opts{s},$opts{a},$aN,$rgb,$scan_out,$track_out,$opts{m});

sub load_chr_size {

    my ($rS,$rH,$file) = @_; 

    print "Loading $file ... ";
    open(my $IN,'<',"$file") or die "$! - Can't open $file\n";
    while(<$IN>){
        chomp;
        my ($chr_name,$chr_size) =  split(/\s+/,$_);

        next if ($chr_size !~ /^\d+$/);

        $rH->{$chr_name} =  $chr_size;
        ${$rS}          += $chr_size;
    }
    close($IN);
    print "done!\n";
}

sub load_sample {

    my ($rS,$rHoHoA,$sample) = @_;

    print "Loading events/$sample ... ";
    my $i = 1;
    open(my $IN,'<',"events/$sample") or die "$! - Can't open events/$sample\n";
    while(<$IN>){
        next if ($_ !~ /^chr/);
        chomp;
        my ($chr,$start,$end,$name,@info)   = split(/\s+/,$_);
        @{$rHoHoA->{$chr}{$i++}}            = ($chr,$start,$end);
    }
    close($IN);
    print "done!\n";

    ${$rS} = $i - 1;
}
sub calc_ss_by_cg{

    my($rH_sbc,$rHoHoA_ebc,$rHoA_aN,$ss_m,$proc,$sample,$adj,$a_by,$scan_out,$genome_size,$events_in_genome) = @_;

    my ($ss_N,$ss_a);
    if ($a_by ne 'chr') {   ($ss_N,$ss_a) = ($genome_size,$events_in_genome);   }

    mkpath("$scan_out/$sample") unless (-e "$scan_out/$sample");

    my $pm = AnyEvent::ForkManager->new(max_workers => $proc);
    $pm->on_start( sub {    my($pm, $pid)           = @_;   } );
    $pm->on_finish(sub {    my($pm, $pid, $status)  = @_;   } );

    print "Scanning chrs ... ";
    open(my $OUT,'>',"$scan_out/$sample/params.txt") or die "$! - Can't open $scan_out/$sample/params.txt\n";
    foreach my $chr ( keys %{$rHoHoA_ebc} ){

        ($ss_N,$ss_a)       = ($rH_sbc->{$chr},scalar(keys %{$rHoHoA_ebc->{$chr}})) if ($a_by eq 'chr');
        @{$rHoA_aN->{$chr}} = ($ss_a,$ss_N);

        print $OUT  "$chr\t$ss_N\t$ss_a\n";

        $pm->start(
                    cb   => \&calc_ss,
                    args => [$ss_N,$ss_a,$ss_m,$chr,$rHoHoA_ebc->{$chr},$sample,$scan_out,$adj]
                );
    }
    close($OUT);

    my $cv = AnyEvent->condvar;
    $pm->wait_all_children(
                            cb => sub {
                                        my ($pm) = @_;
                                        print " done!\n";
                                        $cv->send;
                                      },
                        );

    $cv->recv;
}
sub calc_ss{

    my ($pm,$ss_N,$ss_a,$ss_m,$chr,$rHoA_ebc,$sample,$scan_out,$adj) = @_;

    my (%hash_table,%islands,%info_locus);

    &mk_hash_table(\%hash_table,$ss_a,$ss_m,$ss_N,$ss_m/$ss_N,$chr);
    &look_for_islands($rHoA_ebc,\%islands,\%info_locus,$ss_m,$chr);
    &calc_ss_pvalue(\%islands,\%info_locus,\%hash_table,$ss_m,$chr,$ss_a,$ss_N,$ss_m/$ss_N,$sample,$scan_out,$adj);
    &adjust_pvalue($ss_m,$sample,$chr,$scan_out) if ($adj ne 'un');
    print " $chr";
}

sub mk_hash_table{

    my ($rH_hash_table,$ss_a,$ss_m,$ss_N,$ss_G,$chr) = @_;

    foreach my $ss_k ( 1..$ss_m )   {   $rH_hash_table->{$ss_k} = &ss_pvalue($ss_N,$ss_k,$ss_G,$ss_a,$ss_m);    }
}
sub ss_pvalue{

    my ($ss_N,$ss_k,$ss_G,$ss_a,$ss_m,$chr) = @_;

    my $density_bin  = gsl_ran_binomial_pdf($ss_k,$ss_G,$ss_a);
    return ($ss_k*$ss_N/$ss_m-$ss_a)*$density_bin+2*$density_bin*gsl_sf_hyperg_2F1(1,$ss_k-$ss_a,1+$ss_k,$ss_G/($ss_G-1));
}
sub look_for_islands {
   
    my ($rHoA_ebc,$rHoHoA_i,$rH_il,$ss_m,$chr) = @_;
    
    my ($last_start,$last_end,$island_id,$x);
    
    ($island_id,$x) = (1,1);
    foreach my $i ( sort { $rHoA_ebc->{$a}[1] <=> $rHoA_ebc->{$b}[1] } keys %{$rHoA_ebc} ){

        my ($chr,$start,$end) = @{$rHoA_ebc->{$i}};

        my $info_scan = join("\t",($chr,$start,$end))."\n";

        $rH_il->{$start}++;
    
        if ( $x > 1 ){
            if ( $start - $last_start + 1 > $ss_m ){
                $island_id++;
            }
        }

        push @{$rHoHoA_i->{all}{$island_id}},$info_scan;
        push @{$rHoHoA_i->{pos}{$island_id}},$start;

        $x++;

        $last_start = $start;
    }
}
sub calc_ss_pvalue {

    my ($rHoHoA,$rH_il,$rH_ss,$ss_m,$chr,$ss_a,$ss_N,$ss_G,$sample,$scan_out,$adj) = @_;

    mkpath("$scan_out/$sample/$ss_m/") unless (-e "$scan_out/$sample/$ss_m/");

    my ($OUTL,$OUTS);

    open($OUTL,'>',"$scan_out/$sample/$ss_m/$chr.log") or die "$! - Can't open $scan_out/$sample/$ss_m/$chr.log\n";   
    open($OUTS,'>',"$scan_out/$sample/$ss_m/$chr.txt") or die "$! - Can't open $scan_out/$sample/$ss_m/$chr.txt\n";   

    my ($extend,$status,$middle_point) = (int($ss_m/2),'O',undef);

    foreach my $island_id ( sort {$a <=> $b} keys %{$rHoHoA->{pos}} ){
    
        my @all_pos     = @{$rHoHoA->{pos}{$island_id}};
        my @all_info    = @{$rHoHoA->{all}{$island_id}};

        my ($start_i,$end_i)    = ($all_pos[0],$all_pos[$#all_pos]);
        my ($size)              = ($end_i-$start_i+1);

        if ( $end_i-$start_i < $ss_m ){
            $middle_point       = int(($start_i+$end_i)/2);
            ($start_i,$end_i)   = ($middle_point-$extend,$middle_point+$extend);
            $status = 'U';
        }
        else{
            $status = 'O';
        }
    
        print $OUTL "\n$island_id|$status|$size|$chr:$start_i-$end_i\t" , join("-",@all_pos) , "\n";
        print $OUTL "$island_id|"                                       , join("$island_id|",@all_info);

        my ($start,$end,$c,$inc,$pos);
        for(my $x = $start_i; $x < $end_i; $x += $ss_m){
            $c      = $x + $ss_m;
            $inc    = -1;
            ($start,$end) = ($inc+$x,$inc+$c);
            while( (++$inc < $ss_m) && ($end < $end_i) ){

                ($start,$end)   = ($inc+$x,$inc+$c);
                ($pos)          = ($start);
                print $OUTL "\t$inc:$pos..$end\n";
                my $k;
                while( $pos  <= $end ){
                    if ( exists $rH_il->{$pos} ){
                        $k += $rH_il->{$pos};
                    }
                    $pos++;
                }
                if ( not exists $rH_ss->{$k} ){
                    $rH_ss->{$k} = &ss_pvalue($ss_N,$k,$ss_G,$ss_a,$ss_m);
                }
                print $OUTS "$chr\t$start\t$end\t$k\t$island_id\t$rH_ss->{$k}\n";
            }
        }
    }

    close($OUTS);
    close($OUTL);
}

sub mk_track{

    my ($rA_chr,$rHoA_aN,$sample,$cutoff,$adj,$aN,$rgb,$dir_in,$track_out,$window) = @_;

    mkpath($track_out) unless (-e "$track_out");

    my $hs;
    # open(my $OUT,'>',"$track_out/$sample.$window.$aN.$adj.txt")        or die "$! - Can't open $track_out/$sample.$window.$aN.$adj.txt\n";
    # open(my $COM,'>',"$track_out/$sample.$window.$aN.$adj.graph.txt")  or die "$! - Can't open $track_out/$sample.$window.$aN.$adj.graph.txt\n";
    open(my $HOT,'>',"$track_out/$sample.$window.$aN.$adj.track.bed")  or die "$! - Can't open $track_out/$sample.$window.$aN.$adj.track.bed\n";
    print $HOT "track name=\"$sample\_hsSS.$window\_aN:$aN:$adj\" description=\"$sample\_hs_scanStat_$window\_aN:$aN:$adj\" itemRgb=\"On\"\n";
    
    print "Writing information for hotspots ... ";
    foreach my $chr ( @{$rA_chr} ){

        print " $chr";
        my ($ss_m,$ss_a,$ss_N,$ss_G,$trans_by_chr);

        ($ss_a,$ss_N) = @{$rHoA_aN->{$chr}};

        my $file_in;
        if ($adj eq 'BY')   {   $file_in = "$dir_in/$sample/$window/$adj/$chr.txt";     }
        else                {   $file_in = "$dir_in/$sample/$window/$chr.txt";          }

        if (-s "$file_in"){

            my %island_info;

            &load_island_info(\%island_info,"$dir_in/$sample/$window/$chr.log");

            my (%cluster,%cluster_p);
            &load_scan(\%cluster,\%cluster_p,$file_in,$adj,$cutoff);

            foreach my $island_id ( sort {$a <=> $b} keys %cluster){
                foreach my $cluster_id ( sort {$cluster{$island_id}{$a}[0] <=> $cluster{$island_id}{$b}[0]} keys %{$cluster{$island_id}}){

                    my ($start,$end) = @{$cluster{$island_id}{$cluster_id}}[0..1];
                    my ($c,@breaks) = (0,());
                    foreach my $pos ( @{$island_info{$island_id}} ){
                        if ($pos > $end){   last;   }
                        else{
                            if ($pos >= $start){    push @breaks, $pos;     }
                        }
                        $c++;
                    }

                    my ($s,$e,$k)               = ($breaks[0],$breaks[$#breaks],scalar(@breaks));
                    my ($hot_start,$hot_end,$w) = ($s,$e+1,$e-$s+1);

                    $ss_m = $w;
                    $ss_G = $ss_m/$ss_N;

                    my $pvalue = &ss_pvalue($ss_N,$k,$ss_G,$ss_a,$ss_m);

                    next if ( $pvalue < 0 || $pvalue > $cutoff );

                    my $hs_st_id = "I$island_id\_C$cluster_id";

                    #print $OUT "$hs_st_id\t$chr:$start-$end\t$chr:$hot_start-$hot_end $k $ss_N:$ss_a\n";
                    #print $COM "$chr\t$hot_start\t$hot_end\t$hs_st_id\t$k\t+\t$hot_start\t$hot_end\t$rgb\t$pvalue\n";
                    print $HOT "$chr\t$hot_start\t$hot_end\t$hs_st_id\t$k\t+\t$hot_start\t$hot_end\t$rgb\t$pvalue\n";
                    $hs++;
                }
            }
        }
        # print $OUT "\n";
    }
    # close($OUT);
    # close($COM);
    close($HOT);

    print " done!\n";
    print "Reported $hs hotspots to $track_out/$sample.$window.$aN.$adj.track.bed\n";
}
sub load_island_info{

    my ($rHoA,$file) = @_;

    open(my $LOG,'<',"$file") or die "$! - Can't open $file\n";
    while(<$LOG>){
        if ($_ =~ /^(\d+)\|chr/){
            chomp;
            my @info = split(/\s+/,$_);
            push @{$rHoA->{$1}},$info[1];
        }
    }
    close($LOG);
}
sub load_scan{

    my ($rHoHoA,$rHoHoA_p,$file,$adj,$cutoff) = @_;

    my ($last_start,$last_end,$last_row,$cluster_id);

    $cluster_id = 1;
    open(my $IN,'<',"$file")   or die "$! - Can't open $file\n";
    while(<$IN>){
        chomp;

        if ( $_ =~ /"/ )    {   $_ =~ s/\"//g;  }

        my($tmp,$start,$end,$fr,$island_id,$pvalue) = split(/\s+/,$_);

        next if ( ($adj ne 'BY') && ($pvalue < 0 || $pvalue > $cutoff) );

        if($last_start){
            if( $start >= $last_start && $start < $last_end ){
                $rHoHoA->{$island_id}{$cluster_id}[1]       =  $end;
                push @{$rHoHoA_p->{$island_id}{$cluster_id}},  $pvalue;
            }
            else{
                $cluster_id++;
                @{$rHoHoA->{$island_id}{$cluster_id}} = ($start,$end);
                 push @{$rHoHoA_p->{$island_id}{$cluster_id}},  $pvalue;
            }
        }
        else{
            @{$rHoHoA->{$island_id}{$cluster_id}} = ($start,$end);
            push @{$rHoHoA_p->{$island_id}{$cluster_id}},  $pvalue;
        }
        ($last_start,$last_end,$last_row) = ($start,$end,$_);
    }
    close($IN);
}

sub adjust_pvalue{

    my ($ss_m,$sample,$chr,$scan_out) = @_;

    mkpath("$scan_out/$sample/$ss_m/BY") unless (-e "$scan_out/$sample/$ss_m/BY");

    `Rscript BY.R $ss_m $scan_out $sample $chr`;
}

sub check_args_user{
    
    my ($opts) = @_;

    if ( !$opts->{m} || !$opts->{c} || !$opts->{e} || !$opts->{o} ){
        die "\nHelp:\n",
            "\t-m: window of width m to scan on chromosome;\n",
            "\t-c: chromosome size file;\n",
            "\t-e: events file;\n",
            "\t-o: output directory name;\n",
            "\t-p: max parallel forking count (default: 10);\n",
            "\t-s: significant level (defalut: 0.05);\n",
            "\t-a: adjust p-values using Benjamini & Yekutieli (2001) (default: no).\n\n";
    }
    elsif($opts->{m} && $opts->{m} !~ /^\d+$/){
        die "\n\tERROR: It is expected integer number for 'm' parameter.\n\n";
    }
    elsif($opts->{p} && $opts->{p} !~ /^\d+$/){
        die "\n\tERROR: It is expected integer number for 'p' parameter.\n\n";
    }
    elsif($opts->{s} && ( $opts->{s} < 0 || $opts->{s} > 1 ) ){
        die "\n\tERROR: The 's' parameter must take on values in the interval [0,1].\n\n";
    }
    elsif($opts->{a} && ( $opts->{a} ne 'yes' && $opts->{a} ne 'no') ){
        die "\n\tERROR: It is expected yes or no for 'a' parameter. ($opts->{a})\n\n";
    }
    elsif($opts->{o} && (! -w $opts->{o}) ){
        die "\n\tERROR: '$opts->{o}' directory haven't write permission\n\n";
    }

    $opts->{p} ||= 8;
    $opts->{s} ||= 0.05;
    $opts->{a} ||= 'no';

    if ($opts->{a} eq 'no') {   $opts->{a} = 'un';  }
    else                    {   $opts->{a} = 'BY';  }
}
