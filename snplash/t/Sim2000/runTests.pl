=pod

Run a series of tests on SNPLASH.  The purpose of this suite is to test for changes
between versions rather than to test for actual correctness.

Utilizes a series of inputs called sim2000.{geno, covphen, map}

Compares against a series of outputs called sim2000.canonical.{snpgwa, dandelion, qsnpgwa}

=author Richard T. Guy

=cut


##########
# SNPGWA
##########

system("../src/snpadt -geno sim2000.geno -phen sim2000.covphen -map sim2000.map -out sim2000.test.snpgwa -engine snpgwa ");
$diffs1 = `diff sim2000.test.snpgwa sim2000.canonical.snpgwa`;
$diffs2 = `diff sim2000.test.snpgwa.ref sim2000.canonical.snpgwa.ref`;
$diffs3 = `diff sim2000.test.snpgwa.log sim2000.canonical.snpgwa.log`;

if($diffs1){
	print "Differences exist in SNPGWA output.  To see them, compare sim2000.canonical.snpgwa with sim2000.test.snpgwa\n";
}else{
	print "No differences exist in the SNPGWA output.  \n";
}
