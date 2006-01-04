
use strict;
use warnings;
use Module::Build;
my $build_class = 'Module::Build';
my $custom_build = 'Module::Build::Live';
eval("require $custom_build;");
unless($@) {
	warn "using custom build class:  '$custom_build'\n";
	$build_class = $custom_build;
}

my $builder = $build_class->new(
	module_name         => 'CAD::Calc',
	license             => 'perl',
	dist_version_from   => 'lib/CAD/Calc.pm',
	requires => {
		'Test::More'       => 0,
		'Math::Vec'        => 0.03,
		'Math::Round::Var' => 0.01_01,
		'Math::Complex'    => 0,
		'Math::BigFloat'   => 0,
		#'Math::Geometry::Planar::Offset' => 1.03_01

	},
	add_to_cleanup      => [ qw(CAD-Calc-* META.yml)],
);

$builder->create_build_script();
