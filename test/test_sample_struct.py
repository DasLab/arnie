import sample_structures

sample_seq = 'GGGGAAAACCCC'

def test_sample_seq():
	
	struct_list, ener_list, prob_list = sample_structures.sample_structures(sample_seq, n_samples=10, package='vienna_2')
	print(struct_list)
	print(ener_list)
	print(prob_list)
	return

if __name__=='__main__':
	test_sample_seq()
		#test_pkg_w_bpps(pkg.lower())
