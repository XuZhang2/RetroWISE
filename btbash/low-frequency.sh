python score.py \
	-beam_size 50 \
	-n_best 50 \
	-augmentation 5\
	-targets ./dataset/less-frequency_PtoR_aug5_seed33/test/tgt-test.txt \
	-predictions ./dataset/less-frequency_PtoR_aug5_seed33/test/less-2-tgt.txt \
	-process_number 8\
	-score_alpha 1 \

