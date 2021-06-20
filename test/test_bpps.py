from arnie.bpps import bpps
from arnie.utils import load_package_locations

sample_seq = 'CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG'

def test_bpps(pkg):
	p = bpps(sample_seq, package = pkg)
	print('test bpps %s' % pkg)
	print(p[0])
	return

if __name__=='__main__':
	print("Test: printing first row of bpp matrices")
	package_locs = load_package_locations()
	for pkg in sorted(package_locs.keys()):
		if pkg=='TMP' or pkg.startswith('linear'):
			continue

		test_bpps(pkg.lower())
