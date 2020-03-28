import sys, os, argparse
import arnie.bpps as bpps
from arnie.utils import write_matrix_to_file

if __name__=='__main__':
    p = argparse.ArgumentParser(description=
        """
        Write base pairing probability matrices to files.
        """)
    
    p.add_argument("seq_dir", nargs='+',
                   help="path to dir of *.seq files")
    p.add_argument("-o", help="name of output dir")
    p.add_argument("-p", "--package", default='vienna_2',
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
        bp_matrix = bpps.bpps(seq, package=args.package)
        with open("%s/%s.bpps" % (args.o, seq_id),'w') as f:
            write_matrix_to_file(bp_matrix, f)
