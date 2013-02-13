#!/usr/bin/env perl

use strict;
use warnings;

use PDL;
use PDL::GSL::Bspline;
use PDL::GSL::RNG;

my $rng = PDL::GSL::RNG->new('mt19937');
$rng->set_seed(1);

#my $bspline = PDL::GSL::Bspline::new_bspline(4,10);
#my $bspline_deriv = PDL::GSL::Bspline::new_bspline_deriv(4);

my $N = 500;
my $x = zeroes($N)->xlinvals(0,15.0);
my $y0 = cos($x) * exp(-0.0*$x);
my $sigma = 0.1 * $y0;
my $dy = $rng->ran_gaussian_var($sigma);
my $y = $y0 + $dy;

my $res;

#my ($x1,$y1) = rcols 'splinetest.dat';
my ($x1,$y1) = rcols 'checkcols.dat';

$res = bspline_fit($x1,$y1, { yderiv => 1, ncoeffs => 12, k => 4 });

wcols $x1,$y1,$res->{yfit},$res->{yerr},$res->{yderiv},$res->{yderiv_err}, 'tfit.dat';

#wcols $x,$y0,$y,$res->{yfit},$res->{yerr},
#    $res->{yderiv},$res->{yderiv_err}, $y0-$res->{yfit},"fittest4_12.dat";

wcols $x,$y, 'splinetest.dat';

#$res = bspline_fit($x,$y, { ncoeffs => 17, k => 4 });
#wcols $x,$y0,$y,$res->{yfit},$res->{yerr},"fittest1.dat";

#$res = bspline_fit($x,$y, { yfit => undef, yderiv => 1, ncoeffs => 17, k => 4 });
#wcols $x,$y0,$y,$res->{yderiv},$res->{yderiv_err},"fittest2.dat";
