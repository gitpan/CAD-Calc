package CAD::Calc;
our $VERSION = '0.21';

use Math::Vec qw(NewVec);
# this gets the OffsetPolygon routine (which still needs work)
use Math::Geometry::Planar;
use Math::Geometry::Planar::Offset;
use Math::Complex;
use Math::Round::Var;

use vars qw(
	$linear_precision
	$angular_precision
	$linr
	$angr
	$pi
	);
$linear_precision = 1.0e-7;
$angular_precision = 1.0e-6;
$pi = atan2(1,1) * 4;

require Exporter;
@ISA='Exporter';
@EXPORT_OK = qw (
	pi
	distdivide
	subdivide
	shorten_line
	dist
	dist2d
	line_vec
	slope
	segs_as_transform
	chevron_to_ray
	signdist
	offset
	shift_line
	line_to_rectangle
	isleft
	iswithin
	iswithinc
	unitleft
	unitright
	unit_angle
	angle_reduce
	angle_parse
	angle_quadrant
	collinear
	triangle_angles
	intersection_data
	line_intersection
	seg_line_intersection
	seg_seg_intersection
	seg_ray_intersection
	ray_pgon_int_index
	ray_pgon_closest_index
	perp_through_point
	foot_on_segment
	Determinant
	pgon_as_segs
	pgon_area
	pgon_angles
	pgon_deltas
	pgon_direction
	pgon_bisectors
	sort_pgons_lr
	pgon_splice
	stringify
	pol_to_cart
	cart_to_pol
	print_line
	point_avg
	);


	
use strict;
use Carp;
=pod

=head1 NAME

CAD::Calc - generic cad-related geometry calculations

=head1 AUTHOR

  Eric L. Wilhelm
  ewilhelm at sbcglobal dot net
  http://pages.sbcglobal.net/mycroft

=head1 COPYRIGHT

This module is copyright (C) 2003 by Eric L. Wilhelm and A. Zahner Co.

=head1 LICENSE

This module is distributed under the same terms as Perl.  See the Perl
source package for details.

You may use this software under one of the following licenses:

  (1) GNU General Public License
    (found at http://www.gnu.org/copyleft/gpl.html)
  (2) Artistic License
    (found at http://www.perl.com/pub/language/misc/Artistic.html)

=head1 NO WARRANTY

This software is distributed with ABSOLUTELY NO WARRANTY.  The author
and his employer will in no way be held liable for any loss or damages
resulting from its use.

=head1 Modifications

The source code of this module is made freely available and
distributable under the GPL or Artistic License.  Modifications to and
use of this software must adhere to one of these licenses.  Changes to
the code should be noted as such and this notification (as well as the
above copyright information) must remain intact on all copies of the
code.

Additionally, while the author is actively developing this code,
notification of any intended changes or extensions would be most helpful
in avoiding repeated work for all parties involved.  Please contact the
author with any such development plans.


=cut

=head1 CHANGES

  0.20 
    Added sprintf("%0.9f") to seg_seg_intersection()

  0.21
    Several new functions and features.

=cut
########################################################################

=head1 Configuration

Used to set package global values such as precision.

=cut
########################################################################

=head2 import

Not called directly.  Triggered by the use() function.

  import(%options, @EXPORT_TAGS);

Example:

  use CAD::Calc (
  	-precision => 0.125,
	-angular   => 1.0e-6,
  	qw(
		seg_seg_intersection
		dist2d
		print_line
		)
	);

=cut
sub import {
	## print "import called with @_\n";
	local @ARGV = @_; # shame that Getopt::Long isn't structured better!
	use Getopt::Long;
	GetOptions( '-',
		'precision=f' => \$linear_precision,
		'angular=f'   => \$angular_precision,
		);
	## print "using $linear_precision for linear\n";
	## print "using $angular_precision for angular\n";
	$linr = Math::Round::Var->new($linear_precision);
	$angr = Math::Round::Var->new($angular_precision);
	## print "my linear rounding will be a ", ref($linr), "\n";
	## print "my angular rounding will be a ", ref($angr), "\n";
	CAD::Calc->export_to_level(1, @ARGV);
} # end subroutine import definition
########################################################################

=head1 Constants

=cut
########################################################################

=head2 pi

Returns the value of CAD::Calc::pi

  pi;

=cut
sub pi() {
	return($pi);
} # end subroutine pi definition
########################################################################

=head1 Functions

These are all exported as options.

=cut
########################################################################

=head2 distdivide

Returns a list of point references resulting from dividing $line into
as many parts as possible which are at least $dist apart.

  @points = distdivide(\@line, $dist);

=cut
sub distdivide {
	my($line, $dist) = @_;
	$dist or croak("call to distdivide would cause divide by zero");
	my $ptA = NewVec(@{$line->[0]});
	my $ptB = NewVec(@{$line->[1]});
	my $seg = NewVec($ptB->Minus($ptA));
	my $length = $seg->Length();
	# optionally go for fewer points here?
	my $count = $length / $dist;
	$count = int($count);
	return(subdivide($line, $count));
} # end subroutine distdivide definition
########################################################################

=head2 subdivide

Returns a list of point references resulting from subdividing $line
into $count parts.  The list will be $count-1 items long, (does not
include $line->[0] and $line->[1]);

$line is of the form:  [ [x1, y1, z1], [x2, y2, z2] ] where z1 and z2
are optional.

  @points = subdivide($line, $count);

=cut
sub subdivide {
	my ($line, $count) = @_;
	$count || croak("cannot divide line into zero segments");
	my $ptA = NewVec(@{$line->[0]});
	my $ptB = NewVec(@{$line->[1]});
# 	print "line:  @$ptA -- @$ptB\n";
	my $seg = NewVec($ptB->Minus($ptA));
	my @points;
	for(my $st = 1; $st < $count; $st++) {
		push(@points, [$ptA->Plus( [ $seg->ScalarMult($st / $count) ] ) ] );
		}
	return(@points);
} # end subroutine subdivide definition
########################################################################

=head2 shorten_line

Shortens the line by the distances given in $lead and $tail.

  @line = shorten_line(\@line, $lead, $tail);

=cut
sub shorten_line {
	my ($line, $lead, $tail) = @_;
	my $ptA = NewVec(@{$line->[0]});
	my $ptB = NewVec(@{$line->[1]});
# 	print "line:  @$ptA -- @$ptB\n";
	my $seg = NewVec($ptB->Minus($ptA));
	my $len = $seg->Length();
	($lead + $tail >= $len) && return();
#        croak("CAD::Calc::shorten_line($lead, $tail)\n" .
#                "\t creates inverted line from length: $len\n");
	return(
		[$ptA->Plus([$seg->ScalarMult($lead / $len)])],
		[$ptB->Minus([$seg->ScalarMult($tail / $len)])],
		);
} # end subroutine shorten_line definition
########################################################################

=head2 dist

Returns the direct distance from ptA to ptB.

  dist($ptA, $ptB);

=cut
sub dist {
	my($ptA, $ptB) = @_;
	(ref($ptB) eq "ARRAY") || ($ptB = [0,0,0]);
	my $dist = sqrt(
		($ptB->[0] - $ptA->[0]) ** 2 +
		($ptB->[1] - $ptA->[1]) ** 2 +
		($ptB->[2] - $ptA->[2]) ** 2
		);
	return($dist);
} # end subroutine dist definition
########################################################################

=head2 dist2d

Purposefully ignores a z (2) coordinate.

  dist2d($ptA, $ptB);

=cut
sub dist2d {
	my($ptA, $ptB) = @_;
	# print "ref is: ", ref($ptB), "\n";
	(ref($ptB) eq "ARRAY") || ($ptB = [0,0,0]);
	# print "ptB: @{$ptB}\n";
	my $dist = sqrt(
		($ptB->[0] - $ptA->[0]) ** 2 +
		($ptB->[1] - $ptA->[1]) ** 2
		);
	return($dist);
} # end subroutine dist2d definition
########################################################################

=head2 line_vec

Returns a Math::Vec object representing the vector from $ptA to $ptB
(which is actually a segment.)

  $vec = line_vec($ptA, $ptB);

=cut
sub line_vec {
	return(NewVec(signdist(@_)));
} # end subroutine line_vec definition
########################################################################

=head2 slope

Calculates the 2D slope between points @ptA and @ptB.  Slope is defined
as dy / dx (rise over run.)

If dx is 0, will return the string "inf", which Perl so kindly treats as
you would expect it to (except it doesn't like to answer the question
"what is infinity over infinity?")

  $slope = slope(\@ptA, \@ptB);

=cut
sub slope {
	my @line = @_;
	my @delta = map({$line[1][$_] - $line[0][$_]} 0..1);
	unless($delta[0]) {
		return("inf");
	}
	return($delta[1] / $delta[0]);
} # end subroutine slope definition
########################################################################

=head2 segs_as_transform

Allows two segments to specify transform data.  

Returns: (\@translate, $rotate, $scale), 
	
where:

@translate is a 2D array [$x, $y] basically describing segment @A

$rotate is the angular difference between $A[0]->$B[0] and $A[1]->$B[1]

$scale is the length of $A[1]->$B[1] divided by the length of
$A[0]->$B[0]

  my ($translate, $rotate, $scale) = segs_as_transform(\@A, \@B);

=cut
sub segs_as_transform {
	my ($A, $B) = @_;
	my $av = line_vec(@$A);
	# print_line($A);
	my $sd = line_vec($A->[0], $B->[0]);
	my $ed = line_vec($A->[1], $B->[1]);
	my $ang = $ed->Ang() - $sd->Ang();
	my $sl = $sd->Length();
	$sl or croak("no length for divisor\n");
	my $scale = $ed->Length() / $sl;
	return([$av->[0], $av->[1]], $ang, $scale);
} # end subroutine segs_as_transform definition
########################################################################

=head2 chevron_to_ray

Converts a chevron into a directional line by finding the midpoint
between the midpoints of each edge and connecting to the middle point.

  @line = chevron_to_ray(@pts);

=cut
sub chevron_to_ray {
	my (@pts) = @_;
	(scalar(@pts) == 3) or croak("chevron needs three points");
	my @mids;
	foreach my $seg (0,1) {
		($mids[$seg]) = subdivide([$pts[$seg], $pts[$seg+1]], 2);
	}
	my ($start) = subdivide(\@mids, 2);
	return($start, $pts[1]);
} # end subroutine chevron_to_ray definition
########################################################################

=head2 signdist

Returns the signed distance

  signdist(\@ptA, \@ptB);

=cut
sub signdist {
	my ($ptA, $ptB) = @_;
	my $b = NewVec(@{$ptB});
	return($b->Minus($ptA));
} # end subroutine signdist definition
########################################################################

=head2 offset

Creates a contour representing the offset of @polygon by $dist.
Positive distances are inward when @polygon is ccw.

  @polygons = offset(\@polygon, $dist);

=cut
sub offset {
	my ($polygon, $dist) = @_;
	my @pgons = OffsetPolygon($polygon, $dist);
	return(@pgons);
} # end subroutine offset definition
########################################################################

=head2 intersection_data

Calculates the two numerators and the denominator which are required
for various (seg-seg, line-line, ray-ray, seg-ray, line-ray, line-seg)
intersection calculations.

  ($k, $l, $d) = intersection_data(\@line, \@line);

=cut
sub intersection_data {
	my @l = @_;
	my $n1 = Determinant(
		$l[1][0][0]-$l[0][0][0],
		$l[1][0][0]-$l[1][1][0],
		$l[1][0][1]-$l[0][0][1],
		$l[1][0][1]-$l[1][1][1],
		);
	my $n2 = Determinant(
		$l[0][1][0]-$l[0][0][0],
		$l[1][0][0]-$l[0][0][0],
		$l[0][1][1]-$l[0][0][1],
		$l[1][0][1]-$l[0][0][1],
		);
	my $d  = Determinant(
		$l[0][1][0]-$l[0][0][0],
		$l[1][0][0]-$l[1][1][0],
		$l[0][1][1]-$l[0][0][1],
		$l[1][0][1]-$l[1][1][1],
		);
	return($n1, $n2, $d);

} # end subroutine intersection_data definition
########################################################################

=head2 line_intersection

Returns the intersection point of two lines.

  @pt = line_intersection(\@line, \@line, $tolerance);
  @pt or die "no intersection";

If tolerance is defined, it will be used to sprintf the parallel factor.
Beware of this, it is clunky and might change if I come up with
something better.

=cut
sub line_intersection {
	my @l = (shift, shift);
	my ($tol) = @_;
	foreach my $should (0,1) {
		# print "should have $should\n";
		# print $l[$should], "\n";
		(ref($l[$should]) eq "ARRAY") or warn "not good\n";
	}
	my ($n1, $n2, $d) = intersection_data(@l);
	## print "d: $d\n";
	if(defined($tol)) {
		$d = sprintf("%0.${tol}f", $d);
	}
	if($d == 0) {
		# print "parallel!\n";
		return(); # parallel
	}
	my @pt = (
		$l[0][0][0] + $n1 / $d * ($l[0][1][0] - $l[0][0][0]),
		$l[0][0][1] + $n1 / $d * ($l[0][1][1] - $l[0][0][1]),
		);
#    print "got point: @pt\n";
	return(@pt);
} # end subroutine line_intersection definition
########################################################################

=head2 seg_line_intersection

Finds the intersection of @segment and @line.

  my @pt = seg_line_intersection(\@segment, \@line);
  @pt or die "no intersection";
  unless(defined($pt[1])) {
    die "lines are parallel";
  }

=cut
sub seg_line_intersection {
	my (@l) = @_;
	my ($n1, $n2, $d) = intersection_data(@l);
	# XXX not consistent with line_intersection function
	if(sprintf("%0.9f", $d) == 0) {
		return(0); # lines are parallel
	}
	if( ! (($n1/$d <= 1) && ($n1/$d >=0)) ) {
		return(); # no intersection on segment
	}
	my @pt = (
		$l[0][0][0] + $n1 / $d * ($l[0][1][0] - $l[0][0][0]),
		$l[0][0][1] + $n1 / $d * ($l[0][1][1] - $l[0][0][1]),
		);
	return(@pt);
} # end subroutine seg_line_intersection definition
########################################################################

=head2 seg_seg_intersection

  my @pt = seg_seg_intersection(\@segmenta, \@segmentb);

=cut
sub seg_seg_intersection {
	my (@l) = @_;
	my ($n1, $n2, $d) = intersection_data(@l);
	# print "data $n1, $n2, $d\n";
	if(sprintf("%0.9f", $d) == 0) {
		return(0); # lines are parallel
	}
	if( ! ((sprintf("%0.9f", $n1/$d) <= 1) && (sprintf("%0.9f", $n1/$d) >=0)) ) {
		# warn("n1/d is ", $n1/$d);
		return(); # no intersection on segment a
	}
	if( ! ((sprintf("%0.9f", $n2/$d) <= 1) && (sprintf("%0.9f", $n2/$d) >=0)) ) {
		# warn("n2/d is ", $n2/$d);
		return(); # no intersection on segment b
	}
	my @pt = (
		$l[0][0][0] + $n1 / $d * ($l[0][1][0] - $l[0][0][0]),
		$l[0][0][1] + $n1 / $d * ($l[0][1][1] - $l[0][0][1]),
		);
	return(@pt);
} # end subroutine seg_seg_intersection definition
########################################################################

=head2 seg_ray_intersection

Intersects @seg with @ray, where $ray[1] is the direction of the
infinite ray.

  seg_ray_intersection(\@seg, \@ray);

=cut
sub seg_ray_intersection {
	my (@l) = @_;
	my ($n1, $n2, $d) = intersection_data(@l);
	# $n1 is distance along segment (must be between 0 and 1)
	# $n2 is distance along ray (must be greater than 0)
	if(sprintf("%0.9f", $d) == 0) {
#        print "parallel\n";
		return(0); # lines are parallel
	}
	if( ! (($n1/$d <= 1) && ($n1/$d >=0)) ) {
#        print "nothing on segment\n";
		return(); # no intersection on segment
	}
	if($n2 / $d < 0) {
#        print "nothing on ray\n";
		# segment intersects behind ray
		return();
	}
	my @pt = (
		$l[0][0][0] + $n1 / $d * ($l[0][1][0] - $l[0][0][0]),
		$l[0][0][1] + $n1 / $d * ($l[0][1][1] - $l[0][0][1]),
		);
	return(@pt);

} # end subroutine seg_ray_intersection definition
########################################################################

=head2 ray_pgon_int_index

Returns the first (lowest) index of @polygon which has a segment
intersected by @ray.

  $index = ray_pgon_int_index(\@ray, \@polygon);

=cut
sub ray_pgon_int_index {
	my ($ray, $pgon) = @_;
	(scalar(@$ray) == 2) or croak("not a ray");
	for(my $e = 0; $e < @$pgon; $e++) {
		my $n = $e + 1;
		($n > $#$pgon) && ($n -= @$pgon);
		my $seg = [$pgon->[$e], $pgon->[$n]];
		my @int = seg_ray_intersection($seg, $ray);
		if(defined($int[1])) {
			# print "intersect @int\n";
			return($e);
		}
	}
	return();
} # end subroutine ray_pgon_int_index definition
########################################################################

=head2 ray_pgon_closest_index

Returns the closest (according to dist2d) index of @polygon which has a
segment intersected by @ray.

  $index = ray_pgon_closest_index(\@ray, \@polygon);

=cut
sub ray_pgon_closest_index {
	my ($ray, $pgon) = @_;
	(scalar(@$ray) == 2) or croak("not a ray");
	my @found;
	for(my $e = 0; $e < @$pgon; $e++) {
		my $n = $e + 1;
		($n > $#$pgon) && ($n -= @$pgon);
		my $seg = [$pgon->[$e], $pgon->[$n]];
		my @int = seg_ray_intersection($seg, $ray);
		if(defined($int[1])) {
			# print "intersect @int\n";
			push(@found, [$e, dist2d($ray->[0], \@int)]);
		}
	}
	if(@found) {
		my $least = (sort({$a->[1] <=> $b->[1]} @found))[0];
		return($least->[0]);
	}
	else {
		return();
	}
} # end subroutine ray_pgon_closest_index definition
########################################################################

=head2 perp_through_point

  @line = perp_through_point(\@pt, \@line);

=cut
sub perp_through_point {
	my ($pt, $seg) = @_;
	my @nv = ( # normal vector:
		$seg->[1][1] - $seg->[0][1],
		- ($seg->[1][0] - $seg->[0][0]),
		);
	my @ep = ( # end point of ray
		$pt->[0] + $nv[0],
		$pt->[1] + $nv[1],
		);
	return($pt, \@ep);
} # end subroutine perp_through_point definition
########################################################################

=head2 foot_on_segment

Returns the perpendicular foot of @pt on @seg.  See seg_ray_intersection.

  @pt = foot_on_segment(\@pt, \@seg);

=cut
sub foot_on_segment {
	my ($pt, $seg) = @_;
	return(seg_line_intersection($seg, [perp_through_point($pt, $seg)]));
} # end subroutine foot_on_segment definition
########################################################################

=head2 Determinant

  Determinant($x1, $y1, $x2, $y2);

=cut
sub Determinant {
	my ($x1,$y1,$x2,$y2) = @_;
	return($x1*$y2 - $x2*$y1);
} # end subroutine Determinant definition
########################################################################

=head2 pgon_as_segs

Returns a list of [[@ptA],[@ptB]] segments representing the edges of
@pgon, where segment "0" is from $pgon[0] to $pgon[1]

  @segs = pgon_as_segs(@pgon);

=cut
sub pgon_as_segs {
	my (@pgon) = @_;
	my @segs = ([])x scalar(@pgon);
	for(my $i = -1; $i < $#pgon; $i++) {
		$segs[$i] = [@pgon[$i, $i+1]];
	}
	return(@segs);
} # end subroutine pgon_as_segs definition
########################################################################

=head2 pgon_area

  $area = pgon_area(@polygon);

=cut
sub pgon_area {
	my (@points) = @_;
	my $pgon = Math::Geometry::Planar->new();
	$pgon->points(\@points);
	my $area = $pgon->area();
	return($area);
} # end subroutine pgon_area definition
########################################################################

=head2 pgon_angles

Returns the angle of each edge of polygon in xy plane.  These fall
between -$pi and +$pi due to the fact that it is basically just a call
to the atan2() builtin.

Edges are numbered according to the index of the point which starts
the edge.

  @angles = pgon_angles(@points);

=cut
sub pgon_angles {
	my (@points)  = @_;
	my @angles = (0) x scalar(@points);
	# print "number of angles: @angles\n";
	for(my $i = -1; $i < $#points; $i++) {
		my $vec = NewVec(signdist(@points[$i, $i+1]));
		$angles[$i] = $vec->Ang();
	}
	return(@angles);
} # end subroutine pgon_angles definition
########################################################################

=head2 pgon_deltas

Returns the differences between the angles of each edge of @polygon.
These will be indexed according to the point at which they occur, and
will be positive radians for ccw angles.  Summing the @deltas will yield
+/-2pi (negative for cw polygons.)

  @deltas = pgon_deltas(@pgon);

=cut
sub pgon_deltas {
	my (@pts) = @_;
	my @angles = pgon_angles(@pts);
	return(ang_deltas(@angles));
} # end subroutine pgon_deltas definition
########################################################################

=head2 ang_deltas

Returns the same thing as pgon_deltas, but saves a redundant call to
pgon_angles.

  my @angs = pgon_angles(@pts);
  my @dels = ang_deltas(@angs);

=cut
sub ang_deltas {
	my (@angs) = @_;
	my @deltas = (0)x $#angs;
	for(my $i = 0; $i < @angs; $i++) {
		my $ang = angle_reduce($angs[$i] - $angs[$i-1]);
		$deltas[$i] = $ang;
	}
	return(@deltas);
} # end subroutine ang_deltas definition
########################################################################

=head2 pgon_direction

Returns 1 for counterclockwise and 0 for clockwise.  Uses the sum of the
differences of angles of @polygon.  If this sum is less than 0, the
polygon is clockwise.

  $ang_sum = pgon_direction(@polygon);

=cut
sub pgon_direction {
	my (@pgon) = @_;
	my @angs = pgon_deltas(@pgon);
	return(angs_direction(@angs));
} # end subroutine pgon_direction definition
########################################################################

=head2 angs_direction

Returns the same thing as pgon_direction, but saves a redundant call to
pgon_deltas.

  my @angs = pgon_deltas(@pgon);
  my $dir = angs_direction(@angs);

=cut
sub angs_direction {
	my (@angs) = @_;
	my $sum = 0;
	foreach my $ang (@angs) {
		$sum+= $ang;
	}
	return($sum > 0);
} # end subroutine angs_direction definition
########################################################################

=head2 pgon_bisectors

  pgon_bisectors();

=cut
sub pgon_bisectors {
	warn "unfinished";
	croak "finish it";
} # end subroutine pgon_bisectors definition
########################################################################

=head2 sort_pgons_lr

Sorts polygons by their average points returning a list which reads from
left to right.  (Rather odd place for this?)

  @pgons = sort_pgons_lr(@pgons);

=cut
sub sort_pgons_lr {
	my @pgons = @_;
	# no sense calculating all for naught:
	(scalar(@pgons) > 1) || return(@pgons);
	my @avg;
	foreach my $pgon (@pgons) {
		push(@avg, [point_avg(@$pgon)]);
	}
	my @ord = sort({$avg[$a][0] <=> $avg[$b][0]} 0..$#avg);
	return(@pgons[@ord]);
} # end subroutine sort_pgons_lr definition
########################################################################

=head2 shift_line

Shifts line to right or left by $distance.

  @line = shift_line(\@line, $distance, right|left);

=cut
sub shift_line {
	my ($line, $dist, $dir) = @_;
	my @line = @$line;
	my $mvec;
	if($dir eq "left") {
		$mvec = unitleft(@line);
	}
	elsif($dir eq "right") {
		$mvec = unitright(@line);
	}
	else {
		croak ("direction must be \"left\" or \"right\"\n");
	}
	my @newline = map({[$mvec->Plus($_)]} @line);
	return(@newline);
} # end subroutine shift_line definition
########################################################################

=head2 line_to_rectangle

Creates a rectangle, centered about @line.

  my @rec = line_to_rectangle(\@line, $offset, \%options);

The direction of the returned points will be counter-clockwise around
the original line, with the first point at the 'lower-left' (e.g. if
your line points up, $rec[0] will be below and to the left of
$line[0].)  

=item Available options

  ends => 1|0,   # extend endpoints by $offset (default = 1)

=cut
sub line_to_rectangle {
	my ($ln, $offset, $opts) = @_;
	my %options = (ends => 1);
	(ref($opts) eq "HASH") && (%options = %$opts);
	my @line = @$ln;
	($offset > 0) or
		croak "offset ($offset) must be positive non-zero\n";
	my $a = NewVec(@{$line[0]});
	my $b = NewVec(@{$line[1]});
	# unit vector of line
	my $vec = NewVec(NewVec($b->Minus($a))->UnitVector());
	# crossed with unit vector make unit vector left
	my $perp = NewVec($vec->Cross([0,0,-1]));
	my ($back, $forth);
	if($options{ends}) {
		$back = NewVec($a->Minus([$vec->ScalarMult($offset)]));
		$forth = NewVec($b->Plus([$vec->ScalarMult($offset)]));
	}
	else {
		$back = $a;
		$forth = $b;
	}
	my $left = NewVec($perp->ScalarMult($offset));
	my $right = NewVec($perp->ScalarMult(-$offset));
	# upper and lower here only mean anything
	#  if line originally pointed "up"
	my @ll = $back->Plus($left);
	my @lr = $back->Plus($right);
	my @ur = $forth->Plus($right);
	my @ul = $forth->Plus($left);
	return(\@ll, \@lr, \@ur, \@ul);
} # end subroutine line_to_rectangle definition
########################################################################

=head2 isleft

Returns positive if @point is left of @line.

  isleft(\@line, \@point);

=cut
sub isleft {
	my ($line, $pt) = @_;
	my $isleft = ($line->[1][0] - $line->[0][0]) * 
					($pt->[1] - $line->[0][1]) -
				 ($line->[1][1] - $line->[0][1]) *
				 	($pt->[0] - $line->[0][0]);
	return($isleft);
} # end subroutine isleft definition
########################################################################

=head2 iswithin

Returns true if @pt is within the polygon @bound.

  $fact = iswithin(\@bound, \@pt);

=cut
sub iswithin {
	my ($bnd, $pt) = @_;
	my $winding = 0;
	my @bound = @$bnd;
	for(my $n = -1; $n < $#bound; $n ++) {
		my $next = $n+1;
		my @seg = ($bound[$n], $bound[$next]);
		my $isleft = isleft(\@seg, $pt);
		if($seg[0][1] <= $pt->[1]) {
			if($seg[1][1] > $pt->[1]) {
				($isleft > 0) && $winding++;
#                print "winding up\n";
			}
		}
		elsif($seg[1][1] <= $pt->[1]) {
			($isleft < 0) && $winding--;
#                print "winding up\n";
		}
	} # end for $n
#    print "winding is $winding\n";
	return($winding);
} # end subroutine iswithin definition
########################################################################

=head2 iswithinc

Seems to be consistently much faster than the typical winding-number iswithin.

  iswithinc();

=cut
sub iswithinc {
	my ($bnd, $pt) = @_;
	my $c = 0;
	my @bound = @$bnd;
	my ($x, $y) = @$pt;
	# straight from the comp.graphics.algorithms faq:
	for (my $i = 0, my $j = $#bound; $i < @bound; $j = $i++) {
#        print "checking from $j to $i\n";
		(((($bound[$i][1]<=$y) && ($y<$bound[$j][1])) ||
			(($bound[$j][1]<=$y) && ($y<$bound[$i][1]))) &&
			($x < ($bound[$j][0] - $bound[$i][0]) * ($y - $bound[$i][1]) / ($bound[$j][1] - $bound[$i][1]) + $bound[$i][0]))
			and ($c = !$c);
	}
	return($c);
} # end subroutine iswithinc definition
########################################################################

=head2 unitleft

Returns a unit vector which is perpendicular and to the left of @line.
Purposefully ignores any z-coordinates.

  $vec = unitleft(@line);

=cut
sub unitleft {
	my (@line) = @_;
	my $ln = NewVec(
			NewVec(@{$line[1]}[0,1])->Minus([@{$line[0]}[0,1]])
			);
	$ln = NewVec($ln->UnitVector());
	my $left = NewVec($ln->Cross([0,0,-1]));
	my $isleft = isleft(\@line, [$left->Plus($line[0])]);
# 	print "fact said $isleft\n";
	return($left);
} # end subroutine unitleft definition
########################################################################

=head2 unitright

Negative of unitleft().

  $vec = unitright(@line);

=cut
sub unitright {
	my $vec = unitleft(@_);
	$vec = NewVec($vec->ScalarMult(-1));
	return($vec);
} # end subroutine unitright definition
########################################################################

=head2 unit_angle

Returns a Math::Vec vector which has a length of one at angle $ang (in
the XY plane.) $ang is fed through angle_parse().

  $vec = unit_angle($ang);

=cut
sub unit_angle {
	my ($ang) = @_;
	$ang = angle_parse($ang);
	my $x = cos($ang);
	my $y = sin($ang);
	return(NewVec($x, $y));
} # end subroutine unit_angle definition
########################################################################

=head2 angle_reduce

Reduces $ang (in radians) to be between -pi and +pi.

  $ang = angle_reduce($ang);

=cut
sub angle_reduce {
	my $ang = shift;
	while($ang > $pi) {
		$ang -= 2*$pi;
	}
	while($ang <= -$pi) {
		$ang += 2*$pi;
	}
	return($ang);
} # end subroutine angle_reduce definition
########################################################################

=head2 angle_parse

Parses the variable $ang and returns a variable in radians.  To convert
degrees to radians: $rad = angle_parse($deg . "d")

  $rad = angle_parse($ang);

=cut
sub angle_parse {
	my $ang = shift;
	if($ang =~ s/d$//) {
		$ang *= $pi / 180;
	}
	return($ang);
} # end subroutine angle_parse definition
########################################################################

=head2 angle_quadrant

Returns the index of the quadrant which contains $angle.  $angle is in
radians.

  $q = angle_quadrant($angle);
  @syms = qw(I II III IV);
  print "angle is in quadrant: $syms[$q]\n";

=cut
sub angle_quadrant {
	my $ang = shift;
	my $x = cos($ang);
	my $y = sin($ang);
	my $vert = ($x < 0);
	my $hori = ($y < 0);
	my @list = (
		[0,3],
		[1,2],
		);
	return($list[$vert][$hori]);
} # end subroutine angle_quadrant definition
########################################################################

=head2 collinear

  $fact = collinear(\@pt1, \@pt2, \@pt3);

=cut
sub collinear {
	my @pts = @_;
	(@pts == 3) or croak("must call with 3 points");
	my ($pta, $ptb, $ptc) = @pts;
	my $va = line_vec($pta, $ptb);
	my $vb = line_vec($ptb, $ptc);
	my $cp = NewVec($va->Cross($vb));
	my $ta = $cp->Length();
#    print "my vectors: @$va\n@$vb\n@$cp\n";
#    print "ta: $ta\n";
	return(abs($ta) < 0.001);
} # end subroutine collinear definition
########################################################################

=head2 triangle_angles

Calculates the angles of a triangle based on it's lengths.

  @angles = triangle_angles(@lengths);

The order of the returned angle will be "the angle before the edge".

=cut
sub triangle_angles {
	my @len = @_;
	(@len == 3) or croak("triangle must have 3 sides");
	my @angs = (
		acos(
			($len[2]**2 + $len[0]**2 - $len[1]**2) / 
			(2 * $len[2] * $len[0])
			),
		acos(
			($len[1]**2 + $len[0]**2 - $len[2]**2) / 
			(2 * $len[1] * $len[0])
			),
		);
	$angs[2] = $pi - $angs[0] - $angs[1];
	print "angs: @angs\n";
} # end subroutine triangle_angles definition
########################################################################

=head2 stringify

Turns point into a string rounded according to $rnd.  The optional
$count allows you to specify how many coordinates to use.

  $string = stringify(\@pt, $rnd, $count);

=cut
sub stringify {
	my ($pt, $rnd, $count) = @_;
	# FIXME: # rounding should be able to do fancier things here:
	unless(defined($count)) {
		$count = scalar(@{$pt});
	}
	my $top = $count - 1;
	my $str = join(",", 
		map( {sprintf("%0.${rnd}f", $_)} @{$pt}[0..$top]) );
	return($str);
	
} # end subroutine stringify definition
########################################################################

=head2 pol_to_cart

Convert from polar to cartesian coordinates.

  my ($x, $y, $z) = pol_to_cart($radius, $theta, $z);

=cut
sub pol_to_cart {
	my ($r, $th, $z) = @_;
	my $x = $r * cos($th);
	my $y = $r * sin($th);
	return($x, $y, $z);
} # end subroutine pol_to_cart definition
########################################################################

=head2 cart_to_pol

Convert from polar to cartesian coordinates.

  my ($radius, $theta, $z) = cart_to_pol($x, $y, $z);

=cut
sub cart_to_pol {
	my ($x, $y, $z) = @_;
	my $r = sqrt($x**2 + $y**2);
	my $th = atan2($y, $x);
	return($r, $th, $z);
} # end subroutine cart_to_pol definition
########################################################################

=head2 print_line

  print_line(\@line, $message);

=cut
sub print_line {
	my ($line, $message) = @_;
	unless($message) {
		$message = "line:";
	}
	print join("\n\t", $message, 
		map({join(" ", @$_)} @$line)), "\n";
} # end subroutine print_line definition
########################################################################

=head2 point_avg

Averages the x and y coordinates of a list of points.

	my ($x, $y) = point_avg(@points);

=cut
sub point_avg {
	my(@points) = @_;
	my $i;
	my $num = scalar(@points);
	my $x_avg = 0;
	my $y_avg = 0;
	# print "num is $num\n";
	for($i = 0; $i < $num; $i++) {
		# print "point: $points[$i][0]\n";
		$x_avg += $points[$i][0];
		$y_avg += $points[$i][1];
		}
	# print "avgs:  $x_avg $y_avg\n";
	$x_avg = $x_avg / $num;
	$y_avg = $y_avg / $num;
	return($x_avg, $y_avg); 
} # end subroutine point_avg definition

1;
