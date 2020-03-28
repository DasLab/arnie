import sys, os, argparse
import arnie.bpps as bpps
import numpy as np
from arnie.utils import write_vector_to_file

if __name__=='__main__':
    p = argparse.ArgumentParser(description=
        """Write unpaired posterior probabilities to files.
        """)
    
    p.add_argument("seq_dir", nargs='+',
                   help="path to dir of *.seq files")
    p.add_argument("-o", help="name of output dir")
    p.add_argument("-p", "--package", choices=['vienna_2', 'contrafold_se', 'nupack_95'], default='vienna_2',
                   help="Package to use")

    if len(sys.argv)==1:
        p.print_help(sys.stderr)
        sys.exit(1)

    args = p.parse_args()

    if not os.path.exists('./%s' % args.o):
        os.makedirs('./%s' % args.o)

    for seqfile in args.seq_dir:
        print(seqfile)
        seq=open(seqfile,'r').readlines()[-1].rstrip()
    	seq_id = os.path.basename(seqfile).replace('.seq','')

        unp_vector = 1-np.sum(bpps.bpps(seq, package=args.package),axis=0)

    	with open("%s/%s.unp" % (args.o, seq_id),'w') as f:
    		write_vector_to_file(unp_vector, f)
