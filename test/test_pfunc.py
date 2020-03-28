import pfunc
from utils import load_package_locations_from_yaml
sample_seq = 'GGGGAAAACCCC'

def test_pkg(package):
	
	Z = pfunc.pfunc(sample_seq, package=package, bpps=False)
	print('test %s' % package, Z)
	return

def test_pkg_w_bpps(package):
	
	Z, tmp_file = pfunc.pfunc(sample_seq, package=package, bpps=True)
	print('test %s, tmp file for bpps %s' % (package, tmp_file), Z)
	return

if __name__=='__main__':
	package_locs = load_package_locations_from_yaml('user_default.yaml')
	for pkg in sorted(package_locs.keys()):
		if pkg=='TMP':
			continue
		test_pkg(pkg.lower())
		#test_pkg_w_bpps(pkg.lower())
