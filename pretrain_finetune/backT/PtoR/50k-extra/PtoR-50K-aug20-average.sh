onmt_average_models -output  ./exp/USPTO_extra_50K_PtoR_raw_bt/finetune_average_model_36-40-batch9k-extra-filter5-fp0.55-fullfill.pt \
    -m  exp/USPTO_extra_50K_PtoR_raw_bt/finetune_model.product-reactants-batch9k-extra-filter5-fp0.55-fullfill_step_360000.pt \
        exp/USPTO_extra_50K_PtoR_raw_bt/finetune_model.product-reactants-batch9k-extra-filter5-fp0.55-fullfill_step_370000.pt \
        exp/USPTO_extra_50K_PtoR_raw_bt/finetune_model.product-reactants-batch9k-extra-filter5-fp0.55-fullfill_step_380000.pt \
        exp/USPTO_extra_50K_PtoR_raw_bt/finetune_model.product-reactants-batch9k-extra-filter5-fp0.55-fullfill_step_390000.pt \
        exp/USPTO_extra_50K_PtoR_raw_bt/finetune_model.product-reactants-batch9k-extra-filter5-fp0.55-fullfill_step_400000.pt
