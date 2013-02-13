use PDL::Core::Dev;
use ExtUtils::MakeMaker;
use Data::Dumper;
PDL::Core::Dev->import();

sub get_gsl_libs {
  warn << 'EOW' if ref $PDL::Config{GSL_LIBS};
 The GSL_LIBS config variable must be a string (!)
 not a reference. You should probably leave it undefined
 and rely on gsl-config. Build will likely fail.
EOW

  my $lib = ($PDL::Config{GSL_LIBS} or
	  `gsl-config --libs` or
	  warn "\tno GSL link info (libgsl probably not available)\n");
  my $inc = ($PDL::Config{GSL_INC} or
	  `gsl-config --cflags` or 
	  warn "\tno GSL include info (libgsl probably not available)\n\n");
  chomp $lib; chomp $inc;

  # print STDERR "Lib: $lib\nInc: $inc\n";

  return ($inc,$lib);
}

($GSL_includes, $GSL_libs)  = get_gsl_libs();
   
$package = ["bspline.pd",Bspline,PDL::GSL::Bspline];
%hash = pdlpp_stdargs($package);

$hash{VERSION_FROM} = 'bspline.pd';
$hash{ABSTRACT_FROM} = 'bspline.pd';
$hash{LICENSE} = 'perl';
$hash{AUTHOR} = 'John Lapeyre <jlapeyre@cpan.org>';
push @{$hash{LIBS}} , $GSL_libs;
#push$hash{LIBS} => [ ' -lgsl ' ];
$hash{PREREQ_PM}  = {
    "Exporter" => 0,
    "PDL" => 0
};

sub MY::postamble {pdlpp_postamble($package)};

WriteMakefile( %hash );
