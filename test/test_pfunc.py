from arnie.pfunc import pfunc
from arnie.utils import load_package_locations

sample_seq = 'CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG'

def test_pfunc(package):
	
	Z = pfunc(sample_seq, package=package)
	print('test %s' % package, Z)
	return

if __name__=='__main__':
	package_locs = load_package_locations()
	print(package_locs)
	for pkg in sorted(package_locs.keys()):
		print(pkg)
		if pkg=='TMP':
			continue
		if pkg.startswith('linear'):
			continue

		test_pkg(pkg.lower())