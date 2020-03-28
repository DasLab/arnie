import bpps
import sys
sample_seq = 'GGGGAAAACCCC'

from utils import load_package_locations_from_yaml
def test_bpps(pkg):
	p = bpps.bpps(sample_seq, package = pkg)
	print('test bpps %s' % pkg)
	print(p[0])
	return

if __name__=='__main__':
	print("Test: printing first row of bpp matrices")
	package_locs = load_package_locations_from_yaml('user_default.yaml')
	for pkg in sorted(package_locs.keys()):
		if pkg=='TMP':
			continue

		test_bpps(pkg.lower())
