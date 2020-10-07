#ifndef READFILE_H
#define READFILE_H

/* Raw flight data files from 19.10.2018 */
std::vector<std::string> inputFiles= {
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og/DAMPE_2A_OBS_20181019_20181019T001308_20181019T002610_00000_jEwJChikPnVApmeJU6og.root", 
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T125125_20181019T130429_00000_qEuOl5NzaBT3yrrjYUDk/DAMPE_2A_OBS_20181019_20181019T125125_20181019T130429_00000_qEuOl5NzaBT3yrrjYUDk.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T075207_20181019T080807_00000_uMgc5FmIZB3l69lXjtXj/DAMPE_2A_OBS_20181019_20181019T075207_20181019T080807_00000_uMgc5FmIZB3l69lXjtXj.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T183605_20181019T185324_00000_IJwjQg6lFIW1YmUE3Wyr/DAMPE_2A_OBS_20181019_20181019T183605_20181019T185324_00000_IJwjQg6lFIW1YmUE3Wyr.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T083925_20181019T085324_00000_e0gIu1Qxb3tmE1LHl0ei/DAMPE_2A_OBS_20181019_20181019T083925_20181019T085324_00000_e0gIu1Qxb3tmE1LHl0ei.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T143924_20181019T145940_00000_TMhtd6ymMAIFGOUIqsLI/DAMPE_2A_OBS_20181019_20181019T143924_20181019T145940_00000_TMhtd6ymMAIFGOUIqsLI.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T195721_20181019T201045_00000_laSN2SOCOKgx05ZEvVhG/DAMPE_2A_OBS_20181019_20181019T195721_20181019T201045_00000_laSN2SOCOKgx05ZEvVhG.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T020115_20181019T022028_00000_DgMYdPgNy2V4OswPwHSt/DAMPE_2A_OBS_20181019_20181019T020115_20181019T022028_00000_DgMYdPgNy2V4OswPwHSt.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T230623_20181019T232028_00000_1c9C9bvPsmYVSJQ4Dnhh/DAMPE_2A_OBS_20181019_20181019T230623_20181019T232028_00000_1c9C9bvPsmYVSJQ4Dnhh.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T064635_20181019T070410_00000_Gt4EHFOSqktWxVwz1Eur/DAMPE_2A_OBS_20181019_20181019T064635_20181019T070410_00000_Gt4EHFOSqktWxVwz1Eur.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181020_20181019T235133_20181020T000433_00000_v6lSMJDFPAgkRcM6UDyQ/DAMPE_2A_OBS_20181020_20181019T235133_20181020T000433_00000_v6lSMJDFPAgkRcM6UDyQ.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T211837_20181019T213147_00000_kenVAcYrRCNXQkU3hoi3/DAMPE_2A_OBS_20181019_20181019T211837_20181019T213147_00000_kenVAcYrRCNXQkU3hoi3.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T132547_20181019T134008_00000_uZcuzo59wMtF40R542WL/DAMPE_2A_OBS_20181019_20181019T132547_20181019T134008_00000_uZcuzo59wMtF40R542WL.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T033620_20181019T035502_00000_pswaJ84z4E7kqDAy4y3l/DAMPE_2A_OBS_20181019_20181019T033620_20181019T035502_00000_pswaJ84z4E7kqDAy4y3l.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T100813_20181019T101318_00000_ae1SbppLivUYrbwk8Fp1/DAMPE_2A_OBS_20181019_20181019T100813_20181019T101318_00000_ae1SbppLivUYrbwk8Fp1.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T215525_20181019T220232_00000_xszqYV3l3v38BEnZiAfS/DAMPE_2A_OBS_20181019_20181019T215525_20181019T220232_00000_xszqYV3l3v38BEnZiAfS.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T095540_20181019T100813_00000_WwyCkTl7g41DHK28C1Ru/DAMPE_2A_OBS_20181019_20181019T095540_20181019T100813_00000_WwyCkTl7g41DHK28C1Ru.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T011353_20181019T013310_00000_S6aMpuZxKklgsBundxlh/DAMPE_2A_OBS_20181019_20181019T011353_20181019T013310_00000_S6aMpuZxKklgsBundxlh.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T094233_20181019T095540_00000_vQHVHpeH1XbsznlY29fI/DAMPE_2A_OBS_20181019_20181019T094233_20181019T095540_00000_vQHVHpeH1XbsznlY29fI.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T032316_20181019T033620_00000_jGY6JO15moLJu0nOVCoW/DAMPE_2A_OBS_20181019_20181019T032316_20181019T033620_00000_jGY6JO15moLJu0nOVCoW.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T040954_20181019T042258_00000_d4XNQuIaqIGM0apeyaMy/DAMPE_2A_OBS_20181019_20181019T040954_20181019T042258_00000_d4XNQuIaqIGM0apeyaMy.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T042258_20181019T044218_00000_k7kF0J1cQeDLCg7UhiRw/DAMPE_2A_OBS_20181019_20181019T042258_20181019T044218_00000_k7kF0J1cQeDLCg7UhiRw.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T004635_20181019T010033_00000_Hr8Tq9pFAAV2jkDdX88M/DAMPE_2A_OBS_20181019_20181019T004635_20181019T010033_00000_Hr8Tq9pFAAV2jkDdX88M.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T101359_20181019T102929_00000_Z5dNcBSLi4OkRIumGgp7/DAMPE_2A_OBS_20181019_20181019T101359_20181019T102929_00000_Z5dNcBSLi4OkRIumGgp7.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T114833_20181019T120443_00000_vvrvHSwsl9yd72XF8FEj/DAMPE_2A_OBS_20181019_20181019T114833_20181019T120443_00000_vvrvHSwsl9yd72XF8FEj.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T185405_20181019T190911_00000_u4lyzp7Nv0bDwL71XGH0/DAMPE_2A_OBS_20181019_20181019T185405_20181019T190911_00000_u4lyzp7Nv0bDwL71XGH0.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T030825_20181019T032316_00000_XMreuOv6VYJEWo0ruD62/DAMPE_2A_OBS_20181019_20181019T030825_20181019T032316_00000_XMreuOv6VYJEWo0ruD62.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T121743_20181019T123508_00000_MCQsPuQ8FeSGKo4QxaDA/DAMPE_2A_OBS_20181019_20181019T121743_20181019T123508_00000_MCQsPuQ8FeSGKo4QxaDA.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T150021_20181019T151415_00000_wD9t91C7Av5UX8qx5i2E/DAMPE_2A_OBS_20181019_20181019T150021_20181019T151415_00000_wD9t91C7Av5UX8qx5i2E.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T045829_20181019T051144_00000_cfuc9CfbphpQDOHw3DnE/DAMPE_2A_OBS_20181019_20181019T045829_20181019T051144_00000_cfuc9CfbphpQDOHw3DnE.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T035543_20181019T040954_00000_FPogocJDr4OvecJhgeWu/DAMPE_2A_OBS_20181019_20181019T035543_20181019T040954_00000_FPogocJDr4OvecJhgeWu.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T085324_20181019T090626_00000_0GJ2WYPfnoyWjxCblhUO/DAMPE_2A_OBS_20181019_20181019T085324_20181019T090626_00000_0GJ2WYPfnoyWjxCblhUO.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T173453_20181019T174800_00000_5WjmagLRQmetrzBK9Jox/DAMPE_2A_OBS_20181019_20181019T173453_20181019T174800_00000_5WjmagLRQmetrzBK9Jox.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T194403_20181019T195721_00000_5Iyc1uXSnIyqGeGJguff/DAMPE_2A_OBS_20181019_20181019T194403_20181019T195721_00000_5Iyc1uXSnIyqGeGJguff.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T222954_20181019T225230_00000_XVUkjbWl08yyELyhaQMx/DAMPE_2A_OBS_20181019_20181019T222954_20181019T225230_00000_XVUkjbWl08yyELyhaQMx.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T111658_20181019T113000_00000_2YEKaw1oLIZy0jeJmhQw/DAMPE_2A_OBS_20181019_20181019T111658_20181019T113000_00000_2YEKaw1oLIZy0jeJmhQw.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T053017_20181019T054434_00000_q7GP0XliNDYK6I9BoRwV/DAMPE_2A_OBS_20181019_20181019T053017_20181019T054434_00000_q7GP0XliNDYK6I9BoRwV.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T055735_20181019T061652_00000_YSDHZHRCQG8gxI4w9Dy3/DAMPE_2A_OBS_20181019_20181019T055735_20181019T061652_00000_YSDHZHRCQG8gxI4w9Dy3.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T023520_20181019T024826_00000_VBm2hLFq8GX1I9ROuZ1j/DAMPE_2A_OBS_20181019_20181019T023520_20181019T024826_00000_VBm2hLFq8GX1I9ROuZ1j.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T134008_20181019T135304_00000_fKsUQ5PXAoQUki2iHaxx/DAMPE_2A_OBS_20181019_20181019T134008_20181019T135304_00000_fKsUQ5PXAoQUki2iHaxx.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T071928_20181019T073226_00000_hS8XIcOrPX7hMYW4yvBQ/DAMPE_2A_OBS_20181019_20181019T071928_20181019T073226_00000_hS8XIcOrPX7hMYW4yvBQ.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T120443_20181019T121743_00000_33zgwOjAXhqTiKwKIUA1/DAMPE_2A_OBS_20181019_20181019T120443_20181019T121743_00000_33zgwOjAXhqTiKwKIUA1.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T090626_20181019T092600_00000_smTiDg5FJFpsO3uyCRoI/DAMPE_2A_OBS_20181019_20181019T090626_20181019T092600_00000_smTiDg5FJFpsO3uyCRoI.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T010033_20181019T011353_00000_e4cZGx9YE95aYY74iAZO/DAMPE_2A_OBS_20181019_20181019T010033_20181019T011353_00000_e4cZGx9YE95aYY74iAZO.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T151415_20181019T152717_00000_jCjeqmVJOqJsDfi0cACx/DAMPE_2A_OBS_20181019_20181019T151415_20181019T152717_00000_jCjeqmVJOqJsDfi0cACx.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T044259_20181019T045829_00000_w2NJeaDLLiUmB5IIqRnL/DAMPE_2A_OBS_20181019_20181019T044259_20181019T045829_00000_w2NJeaDLLiUmB5IIqRnL.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T061733_20181019T063321_00000_AnX2RDBRkNE4rk6HIx4k/DAMPE_2A_OBS_20181019_20181019T061733_20181019T063321_00000_AnX2RDBRkNE4rk6HIx4k.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T204310_20181019T205616_00000_C3PnXCV7GEG7QNNSnkOi/DAMPE_2A_OBS_20181019_20181019T204310_20181019T205616_00000_C3PnXCV7GEG7QNNSnkOi.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T080807_20181019T082116_00000_oC6znBnNRu2bPK9HEonF/DAMPE_2A_OBS_20181019_20181019T080807_20181019T082116_00000_oC6znBnNRu2bPK9HEonF.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T232028_20181019T233706_00000_cu6fCNZZzBuq1VAgiuIw/DAMPE_2A_OBS_20181019_20181019T232028_20181019T233706_00000_cu6fCNZZzBuq1VAgiuIw.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T174800_20181019T180848_00000_oF2kmFs8kRNQslvtb3WU/DAMPE_2A_OBS_20181019_20181019T174800_20181019T180848_00000_oF2kmFs8kRNQslvtb3WU.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T214535_20181019T215525_00000_AZaiD5sQtFVnQoHOUmSc/DAMPE_2A_OBS_20181019_20181019T214535_20181019T215525_00000_AZaiD5sQtFVnQoHOUmSc.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T073226_20181019T075126_00000_j4z4Kjo216P6eX9QAcau/DAMPE_2A_OBS_20181019_20181019T073226_20181019T075126_00000_j4z4Kjo216P6eX9QAcau.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T013351_20181019T014809_00000_XIFZ9eLX2zKhysA3L1kw/DAMPE_2A_OBS_20181019_20181019T013351_20181019T014809_00000_XIFZ9eLX2zKhysA3L1kw.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T225311_20181019T230623_00000_lFk0eNCeKEu1zu948J22/DAMPE_2A_OBS_20181019_20181019T225311_20181019T230623_00000_lFk0eNCeKEu1zu948J22.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T110115_20181019T111658_00000_RckbOL5uXuW2WfCFqRRx/DAMPE_2A_OBS_20181019_20181019T110115_20181019T111658_00000_RckbOL5uXuW2WfCFqRRx.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T192220_20181019T194321_00000_YMqzVuV8IJ193dbEKfC6/DAMPE_2A_OBS_20181019_20181019T192220_20181019T194321_00000_YMqzVuV8IJ193dbEKfC6.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T014809_20181019T020115_00000_M5CH3a0E9Vsi6XzZcNyk/DAMPE_2A_OBS_20181019_20181019T014809_20181019T020115_00000_M5CH3a0E9Vsi6XzZcNyk.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T163455_20181019T164823_00000_1UWTHBKIC3bWIujsTk0f/DAMPE_2A_OBS_20181019_20181019T163455_20181019T164823_00000_1UWTHBKIC3bWIujsTk0f.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T131657_20181019T132506_00000_JQjozGwrEZ39TH3NtPEc/DAMPE_2A_OBS_20181019_20181019T131657_20181019T132506_00000_JQjozGwrEZ39TH3NtPEc.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T130429_20181019T131656_00000_Taeg2nNNCBdSmnBUdVVb/DAMPE_2A_OBS_20181019_20181019T130429_20181019T131656_00000_Taeg2nNNCBdSmnBUdVVb.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T205616_20181019T211756_00000_EE0eGCbDsrkgg1euBXNW/DAMPE_2A_OBS_20181019_20181019T205616_20181019T211756_00000_EE0eGCbDsrkgg1euBXNW.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T063321_20181019T064635_00000_o70pm3N5GRY0Luvq2yhw/DAMPE_2A_OBS_20181019_20181019T063321_20181019T064635_00000_o70pm3N5GRY0Luvq2yhw.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T160027_20181019T161333_00000_ldG9PQoWHTQV7kat40YU/DAMPE_2A_OBS_20181019_20181019T160027_20181019T161333_00000_ldG9PQoWHTQV7kat40YU.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T180929_20181019T182249_00000_t9WshIcZxBgrlXEEETC8/DAMPE_2A_OBS_20181019_20181019T180929_20181019T182249_00000_t9WshIcZxBgrlXEEETC8.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T182249_20181019T183605_00000_Th1YVRCVUOOl1IfNt9c6/DAMPE_2A_OBS_20181019_20181019T182249_20181019T183605_00000_Th1YVRCVUOOl1IfNt9c6.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T102929_20181019T104233_00000_FrfWq3FSRSDzEeHkdCux/DAMPE_2A_OBS_20181019_20181019T102929_20181019T104233_00000_FrfWq3FSRSDzEeHkdCux.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T190911_20181019T192220_00000_hQ96i0ZPZaajAOsOVi3N/DAMPE_2A_OBS_20181019_20181019T190911_20181019T192220_00000_hQ96i0ZPZaajAOsOVi3N.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T202136_20181019T202758_00000_lRSy55QELnA8uA5tquK0/DAMPE_2A_OBS_20181019_20181019T202136_20181019T202758_00000_lRSy55QELnA8uA5tquK0.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T164823_20181019T170130_00000_A4QnM7ssOB1yoy8IxUpf/DAMPE_2A_OBS_20181019_20181019T164823_20181019T170130_00000_A4QnM7ssOB1yoy8IxUpf.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T082116_20181019T083448_00000_P56Pf3lDn0p67IM5htGm/DAMPE_2A_OBS_20181019_20181019T082116_20181019T083448_00000_P56Pf3lDn0p67IM5htGm.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T202839_20181019T204310_00000_WvrKD8hVG9Jpyg484hxa/DAMPE_2A_OBS_20181019_20181019T202839_20181019T204310_00000_WvrKD8hVG9Jpyg484hxa.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T083448_20181019T083844_00000_8eThNYJz29Ur02yk65uI/DAMPE_2A_OBS_20181019_20181019T083448_20181019T083844_00000_8eThNYJz29Ur02yk65uI.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T171931_20181019T173453_00000_25bMENVZFP77MFDKO0AC/DAMPE_2A_OBS_20181019_20181019T171931_20181019T173453_00000_25bMENVZFP77MFDKO0AC.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T142623_20181019T143924_00000_7iYL9YNfQ73PfZoQTlg9/DAMPE_2A_OBS_20181019_20181019T142623_20181019T143924_00000_7iYL9YNfQ73PfZoQTlg9.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T022109_20181019T023520_00000_DaNAFi3Z8GKcGjDw5uDa/DAMPE_2A_OBS_20181019_20181019T022109_20181019T023520_00000_DaNAFi3Z8GKcGjDw5uDa.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T141123_20181019T142623_00000_MV71zT0sryHPnddScFAC/DAMPE_2A_OBS_20181019_20181019T141123_20181019T142623_00000_MV71zT0sryHPnddScFAC.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T092641_20181019T094233_00000_6LIo6dNat1E22cdKmo0q/DAMPE_2A_OBS_20181019_20181019T092641_20181019T094233_00000_6LIo6dNat1E22cdKmo0q.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181020_20181019T233747_20181019T235133_00000_MoJCW8innKJKevlfrOE4/DAMPE_2A_OBS_20181020_20181019T233747_20181019T235133_00000_MoJCW8innKJKevlfrOE4.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T114137_20181019T114752_00000_4txRWjHxoVirZMxp0Vwb/DAMPE_2A_OBS_20181019_20181019T114137_20181019T114752_00000_4txRWjHxoVirZMxp0Vwb.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T113000_20181019T114137_00000_l0W187OCMMz3AnQuFnGi/DAMPE_2A_OBS_20181019_20181019T113000_20181019T114137_00000_l0W187OCMMz3AnQuFnGi.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T002610_20181019T004554_00000_OWsUvv6q4hnMtBE09fYu/DAMPE_2A_OBS_20181019_20181019T002610_20181019T004554_00000_OWsUvv6q4hnMtBE09fYu.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T070451_20181019T071928_00000_Ikj9nbuXbEQyVtdF1k7n/DAMPE_2A_OBS_20181019_20181019T070451_20181019T071928_00000_Ikj9nbuXbEQyVtdF1k7n.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T051144_20181019T052936_00000_7gZ7kxWoP7dGP8Yx181X/DAMPE_2A_OBS_20181019_20181019T051144_20181019T052936_00000_7gZ7kxWoP7dGP8Yx181X.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T152717_20181019T154416_00000_eH2bKfbqia65LC3GXbb4/DAMPE_2A_OBS_20181019_20181019T152717_20181019T154416_00000_eH2bKfbqia65LC3GXbb4.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T054434_20181019T055735_00000_o0sMDCt3NoHxFfhDOmey/DAMPE_2A_OBS_20181019_20181019T054434_20181019T055735_00000_o0sMDCt3NoHxFfhDOmey.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T161333_20181019T163414_00000_uYsx85EYSzyvAtVWhJUm/DAMPE_2A_OBS_20181019_20181019T161333_20181019T163414_00000_uYsx85EYSzyvAtVWhJUm.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T201046_20181019T202136_00000_5Xdeo3VwP6EAPjusasbA/DAMPE_2A_OBS_20181019_20181019T201046_20181019T202136_00000_5Xdeo3VwP6EAPjusasbA.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T123549_20181019T125125_00000_at1hSoidQ4qRo3gnRivT/DAMPE_2A_OBS_20181019_20181019T123549_20181019T125125_00000_at1hSoidQ4qRo3gnRivT.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T170130_20181019T171850_00000_gkstdBloKQgPsbKq7fqu/DAMPE_2A_OBS_20181019_20181019T170130_20181019T171850_00000_gkstdBloKQgPsbKq7fqu.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T154457_20181019T160027_00000_hNy0Rw59lbjJM4kIftdu/DAMPE_2A_OBS_20181019_20181019T154457_20181019T160027_00000_hNy0Rw59lbjJM4kIftdu.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T220313_20181019T221659_00000_HDQepXyrFOIFd41ayivB/DAMPE_2A_OBS_20181019_20181019T220313_20181019T221659_00000_HDQepXyrFOIFd41ayivB.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T221659_20181019T222954_00000_MDDrumSKAvusrM4VpiEa/DAMPE_2A_OBS_20181019_20181019T221659_20181019T222954_00000_MDDrumSKAvusrM4VpiEa.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T135304_20181019T140942_00000_nEqbocZlrLBnmBbKLvkr/DAMPE_2A_OBS_20181019_20181019T135304_20181019T140942_00000_nEqbocZlrLBnmBbKLvkr.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T024826_20181019T030744_00000_klN4a8WoJqTCLnJib1Yw/DAMPE_2A_OBS_20181019_20181019T024826_20181019T030744_00000_klN4a8WoJqTCLnJib1Yw.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T213147_20181019T214535_00000_x7fx4vPZFCX0ahpGVVyk/DAMPE_2A_OBS_20181019_20181019T213147_20181019T214535_00000_x7fx4vPZFCX0ahpGVVyk.root",
"/beegfs/dampe/prod/FM/FlightData/2A/20181019/DAMPE_2A_OBS_20181019_20181019T104233_20181019T110034_00000_Sv0xXZyhDzXsUV55gMie/DAMPE_2A_OBS_20181019_20181019T104233_20181019T110034_00000_Sv0xXZyhDzXsUV55gMie.root"
};

#endif // READFILE_H

