install:
	formatdb -i trans_peps.fa -p T
	hmmpress trans_models.hmm
