import os
from arnie.utils import load_package_locations


def test_settings():
    package_locs = load_package_locations()
    for k in package_locs.keys():
        print(k)
        assert os.path.isdir(package_locs[k])
    return


if __name__ == '__main__':
    test_settings()
