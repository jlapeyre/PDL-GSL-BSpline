use PDL::Core::Dev;
use ExtUtils::MakeMaker;
use Data::Dumper;
PDL::Core::Dev->import();
   
$package = ["bspline.pd",Bspline,PDL::GSL::Bspline];
%hash = pdlpp_stdargs($package);

$hash{VERSION_FROM} = 'bspline.pd';
$hash{ABSTRACT_FROM} = 'bspline.pd';
$hash{LICENSE} = 'perl';
$hash{AUTHOR} = 'John Lapeyre <jlapeyre@cpan.org>';
push @{$hash{LIBS}} , ' -lgsl ';
#push$hash{LIBS} => [ ' -lgsl ' ];
$hash{PREREQ_PM}  = {
    "Exporter" => 0,
    "PDL" => 0
};

sub MY::postamble {pdlpp_postamble($package)};

WriteMakefile( %hash );
