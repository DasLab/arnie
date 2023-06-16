from arnie.pfunc import pfunc
from arnie.utils import load_package_locations

sample_seq = 'CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG'


def test_pfunc(package):

    Z = pfunc(sample_seq, package=package)
    print('test %s' % package, Z)
    return


if __name__ == '__main__':
    package_locs = load_package_locations()
    for pkg in sorted(package_locs.keys()):

        if (pkg == 'TMP') or (
            pkg.startswith('linear')) or (
            pkg in ['hotknots', 'ipknot', 'knotty', 'pknots', 'spotrna',
                    'spotrna_conda_env', 'e2efold', 'e2efold_conda_env',
                    'spotrna2']):
            print(f'{pkg} not tested.')
            continue
        print(pkg)
        test_pfunc(pkg.lower())

