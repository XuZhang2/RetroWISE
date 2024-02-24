onmt_train -config pretrain_finetune/backT/PtoR/50k-extra/PtoR-50K-aug20-config.yml
bash pretrain_finetune/backT/PtoR/50k-extra/PtoR-50K-aug20-average.sh
onmt_translate -config pretrain_finetune/backT/PtoR/50k-extra/PtoR-50K-aug20-translate.yml
bash btbash/score-extra-50k.sh
