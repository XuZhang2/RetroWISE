python score.py \
	-beam_size 50 \
	-n_best 50 \
	-augmentation 20 \
	-predictions ./exp/USPTO_extra_50K_PtoR_raw_bt/finetune_average_model_36-40-results-batch9k-extra-filter5-fp0.55-fullfill.txt \
	-targets ./dataset/USPTO_50K_PtoR_aug20/test/tgt-test.txt \
	-process_number 64 \
	-score_alpha 1 \
