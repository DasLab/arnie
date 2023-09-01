from arnie import pfunc

# 3MXH c-di-GMP riboswitch, has coaxial stacking
sample_seq = 'GGUCACGCACAGGGCAAACCAUUCGAAAGAGUGGGACGCAAAGCCUCCGGCCUAAACCAUUGCACUCCGGUAGGUAGCGGGGUUACCGAUGG'


def test_pkg(package, coaxial=True):

    Z = pfunc.pfunc(sample_seq, package=package, bpps=False, coaxial=coaxial)
    print('test %s' % package, Z)
    return None


if __name__ == '__main__':
    for pkg in ['vfold_0', 'vfold_1']:
        for coaxial in [True, False]:
            print(pkg, "coaxial %d" % coaxial)
            test_pkg(pkg, coaxial=coaxial)
