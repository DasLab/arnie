from arnie.free_energy import free_energy
from arnie.mfe import mfe

seq = 'CGCUGUCUGUACUUGUAUCAGUACACUGACGAGUCCCUAAAGGACGAAACAGCG'
dG = free_energy(seq, linear=True, DEBUG=True)
print(dG)

dG = free_energy(seq, linear=True, package='contrafold', DEBUG=True)
print(dG)

dG = free_energy(seq, linear=True, package='eternafold', DEBUG=True)
print(dG)

struct = mfe(seq, linear=True)
print(struct)
struct = mfe(seq, linear=True, package='contrafold')
print(struct)
struct = mfe(seq, linear=True, package='eternafold')
print(struct)

