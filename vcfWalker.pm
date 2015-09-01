#!/usr/bin/perl
#GJR 9/1/2015
# vcfWalker.pm
# a package providing various parsing and printing functions for vcf files
# current methods are 'new' (constructor) and 'walk'
# internal subs include 'printBimbam'



package vcfWalker;    # This is the Class;
use strict;
use Scalar::Util 'blessed';


#constructor for vcfWalker
#possible paramaters:
#  -vcfFile, -sampleFile, -locusFile, -fields

sub new {
    my ($caller, @args) = @_;
    my $class = ref($caller) || $caller;
    

    my %params = @args;
    @params{ map { lc $_ } keys %params } = values %params; # lowercase keys

    unless( defined $params{-vcfFile} ) {
      die "Must specify -vcfFile (OPTIONAL: -sampleFile,-locusFile)";
    }
    my $self = {};
    bless($self,$class);
    if (defined $params{-locusFile} ) {
      $self->{locusFile} = $params{-locusFile};
      my %locusHash;
      open IN, "<$self->{locusFile}";
      while (<IN>) {
	chomp;
	$locusHash{$_}=0;
      }
      close IN;
      $self->{locusHashRef}=\%locusHash;
    }
    if (defined $params{-sampleFile} ) {
      $self->{sampleFile} = $params{-sampleFile};
      my %idhash;
      open IN, "<$self->{sampleFile}";
      while (<IN>) {
	chomp;
	my @vals=split "\t";
	if (@vals > 1) {$idhash{$vals[0]}=$vals[1]} else {$idhash{$vals[0]} = 0} ;
      }
      close IN;
      $self->{sampleHashRef}=\%idhash;
    }
      
    $self->{vcfFile} = $params{-vcfFile};
    if (defined $params{-fields} ) { $self->{fieldsRef} = \$params{-fields} };
    if (defined $params{-outfile} ) { $self->{outfile} = $params{-outfile} };
    return $self;
 
}


sub walk {
  my $self=shift;
  die "not a vcfWalker object" unless blessed $self == 'vcfWalker';
  my %params = @_;
  @params{ map { lc $_ } keys %params } = values %params; # lowercase keys

  if (defined $params{-function} ) {$self->{function} = $params{-function}}
  # -print can may be 'vcf' or 'delim'
  if (defined $params{-print} ) {$self->{print} = $params{-print}} else {$self->{print} = 'no'}
  if (defined $params{-store} ) {$self->{store} = $params{-store}} else {$self->{store} = 'no'}
  if (defined $params{-header} ) {$self->{header} = $params{-header}} else {$self->{header} = 'false'}
  if (defined $params{-delim} ) {$self->{delim} = $params{-delim}} else {$self->{delim} = "\t"}

  my $d = $self->{delim};
  $"=$d;
  #walk through vcf file
  open IN, "<$self->{vcfFile}";
  open OUT, ">$self->{outfile}" if $params{print} !~ m/no/;
  print OUT "$self->{fields}\n" unless $self->{header} =~ m/false/i;
  my @retainInd;
  my @ids;
  my $readInfo="F";
  my %fieldsInds;
  while (<IN>) {

    ######## process vcf headers ###########
    if (/^\#\#/) {
      print OUT "$_" if $self->{print} =~ m/vcf/; #print 'vcf'
      next;
    }
    chomp;
    my @vals = split "\t";

    ######## process column vcf names ###########
    if (/^\#/) {
      if ($self->{print} =~ m/vcf/) { #print 'vcf'
	local $"="\t";
	print OUT "@vals[0..8]";
      }
      if (defined $self->{sampleHashRef}) {
	@retainInd=();
	@ids=();
	my $ind=9;
	for my $id (@vals[9..$#vals]) {
	  if (exists $self->{sampleHashRef}->{$id}) {
	    push @retainInd, $ind;
	    push @ids, $id;
	    print OUT "\t$id" unless $self->{print} =~ m/no/; #print 'vcf' or 'delim' 
	  } 
	  $ind++;
	}
	print OUT "\n";
      } else {
	@retainInd=(9..$#vals);
	@ids=@vals[9..$#vals];
	local $"="\t";
	print OUT "\t@vals[9..$#vals]\n" if $self->{print} =~ m/vcf/;
      }
      next;
    }


    ######## process vcf data ###########
    my $locus="$vals[0]_$vals[1]";
    #skip locus if locus file provided and provided file does not contain the locus
    next if defined $self->{locusHashRef} and not exists $self->{locusHashRef}->{$locus};

    #print first 9 descriptor columns if printing to vcf
    if ($self->{print} =~ m/vcf/) { #print 'vcf'
      local $"="\t";
      print OUT "@vals[0..8]";
    }


    #print all if print = 'vcf' and there is no sample file provided
    if ($self->{print} =~ m/vcf/ and not defined $self->{sampleHashRef} ) {
      local $"="\t";
      print OUT "\t@vals[9..$#vals]\n";
      next;
    }

    #if first line of data, capture the fields and indices in %fieldsInds
    if ($readInfo =~ m/F/i ) {
      my $ind=0;
      my @fields = split "\:", @vals[8];
      for my $field (@fields) {
	$fieldsInds{$field} = $ind;
	$ind++;
      }
      $readInfo="T";
    }

    #process data columns for bimbam format if print = 'bimbam' selected
    #for now, print = 'bimbam' assumes that a 'PL' field is present
    if ($self->{print} =~ m/bimbam/) {
      printBimbam($_,$fieldsInds{'PL'},@retainInd);
      next;
    }
    
    #process data columns for 'delim' or 'vcf' print format
    #selctivity for sample output determined by @retainInd
    for my $info (@vals[@retainInd]) {

      #if missing values for genotype information (data), possibly print and skip subsequent steps
      if ($info =~ m/\.\/\./) {
	print OUT "\t$info" if $self->{print} =~ m/vcf/;
	print OUT "\tNA" if $self->{print} =~ m/delim/;
	next;
      }

      #if print = 'vcf', print and skip subsequent steps
      if ($self->{print} =~ m/vcf/) {
	print OUT "\t$info";
	next;
      }

      #process data if selecting fields for delimited printing
      my @AllInfo=split ":", $info;
      my @fieldVals;
      #process for 'PL' field
      if ( my $ind = $fieldsInds{'PL'} ) {
	my @fieldVals=split ",",$AllInfo[$ind];
	#convert from phred scale to probability
	@fieldVals = map {10**(-$_/10)} @fieldVals if exists $self->{function} and $self->{function} =~ m/phredToLik/;
	#re-scale based on C*(p1+p2+p3)=1 (i.e., scaled probabilities sum to one)
	@fieldVals = map {1/(eval join '+', @fieldVals)*$_} @fieldVals if exists $self->{function} and $self->{function} =~ m/likToProb/;
	print OUT "\t@fieldVals";
      }
    }
    print OUT "\n";
 
    



    
  }
  close OUT if $params{print} !~ m/no/;
  close IN;
}





##### non-object subs ###################

sub printBimbam {
  my $line=shift;
  my $fieldInd-shift;
  my @vals = split "\t", $line;
  my $locus="$vals[0]_$vals[1]";
  my @meanGenos;
  my @nonMissing;
  for my $info (@_) {
    if ($info =~ m/\.\/\./) {
      push @meanGenos, "NA";
      next;
    }
    my @AllInfo=split ":", $info;
    my @PLs=split ",",$AllInfo[$fieldInd];
    #convert from phred scale to probability
    @PLs = map {10**(-$_/10)} @PLs;
    #re-scale based on C*(p1+p2+p3)=1 (i.e., scaled probabilities sum to one)
    @PLs = map {1/(eval join '+', @PLs)*$_} @PLs;
    my $meanGeno = (2*$PLs[2]+$PLs[1]);
    push @meanGenos, $meanGeno;
    push @nonMissing, $meanGeno;
  }
  #hard coded threshold of max 40% missing values
  next if (@nonMissing/@meanGenos) < 0.4;
  print OUT "$locus,$vals[4],$vals[3],";
  my $mean = (eval join '+', @nonMissing)/@nonMissing;
  @meanGenos = replace($mean,@meanGenos);
  print OUT "@meanGenos\n";
}



chdir '/homes/gragland/shared/Rad_data/Rpom/ScottSelection';
my @fields=('PL');

my $newobj=vcfWalker->new(-vcfFile => 'snps.ScottSelection.05-400-20.small.vcf',
			  #-fields => \@fields,
			  -locusFile => 'testsites.txt',
			  -sampleFile => 'sample.file',
			  -outfile => 'out.vcf'
			 );

#my @functions=('phredToLik',likToProb)

my $var = $newobj->walk(#-function => \@functions,
			-print => 'vcf',
			-store => 'no'
			
		       );



my $dum;
