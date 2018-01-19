function exp_dict = T1_datalist()

exp_dict = containers.Map;
%exp_dict('2015_02_19') = {'adch_62f','adch_62g','adch_62l','adch_63h','adch_71b','adch_73h','adch_82f','adch_82k','adch_83d','adch_83e','adch_83f','adch_83g','adch_83h'};
%exp_dict('2015_02_26') = {'adch_24f','adch_25d','adch_34h','adch_36a','adch_36b','adch_46a','adch_46i','adch_46k'};
%exp_dict('2015_03_11') = {'adch_62j','adch_82a','adch_82b','adch_82h','adch_83i','adch_83j'};

exp_dict('2015_06_18') = {'adch_33a','adch_33b','adch_44f'};
%exp_dict('2015_07_28') = {'adch_22b','adch_22c','adch_41e'};

%exp_dict('2016_09_08_L1') = {'adch_75d'}; 
%exp_dict('2016_09_08_R1') = {'adch_21e','adch_22b'};

%% NOT WORKING.

%exp_dict('2017_03_07_L1') = {'adch_34c','adch_46b'}; % NOT WORKING. No STA is extracted
%exp_dict('2017_03_07_R1') = {'adch_48g'};% NOT WORKING. No STA is extracted

%exp_dict('2017_03_16_R1') = {'adch_36b'};% NOT WORKING. It has 36 trials but 54 rexp files

%exp_dict('2017_03_21') = {'adch_56g','adch_67e'};% NOT WORKING. It has 36 trials but 54 rexp files
%exp_dict('2017_04_12') = {'adch_52b'};% NOT WORKING. That cell is not existing
%exp_dict('2017_04_14') = {'adch_23d','adch_33g','adch_44e'};

end