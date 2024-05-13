import argparse
import json
import re

parser = argparse.ArgumentParser(description='Extract integrals into json shape variations.')
parser.add_argument('--input', required=True, type=str, help="input root file")
parser.add_argument('--output', required=True, type=str, help="output json file")
parser.add_argument('--year', required=False, default=2016, type=int, help="bin")
args = parser.parse_args()

processes = {
	 "DY" : "DY",
	 #"DY_LM" : "DY_LM",
	 "EWK" : "EWK",
	 "TT" : "TT",
	 #"TTX" : "TTX",
	 "TW" : "TW",
	 "VV" : "VV",
	 "VVV" : "VVV",
	 "W" : "W",
	 "WH_htt" : "VH_htt",
	 "ZH_hbb" : "VH_hbb",
	 "ggHH_kl_1_kt_1_hbbhtt" : "ggHH_kl_1_kt_1_hbbhtautau",
	 "ggHH_kl_2p45_kt_1_hbbhtt" : "ggHH_kl_2p45_kt_1_hbbhtautau",
	 "ggHH_kl_5_kt_1_hbbhtt" : "ggHH_kl_5_kt_1_hbbhtautau",
	 "ggH_htt" : "ggH_htt",
	 #"qqHH_CV_1_C2V_0_kl_1_hbbhtt" : "qqHH_CV_1_C2V_0_kl_1_hbbhtautau",
	 "qqHH_CV_1_C2V_1_kl_0_hbbhtt" : "qqHH_CV_1_C2V_1_kl_0_hbbhtautau",
	 "qqHH_CV_1_C2V_1_kl_1_hbbhtt" : "qqHH_CV_1_C2V_1_kl_1_hbbhtautau",
	 "qqHH_CV_1_C2V_1_kl_2_hbbhtt" : "qqHH_CV_1_C2V_1_kl_2_hbbhtautau",
	 "qqHH_CV_1_C2V_2_kl_1_hbbhtt" : "qqHH_CV_1_C2V_2_kl_1_hbbhtautau",
	 "qqHH_CV_1p5_C2V_1_kl_1_hbbhtt" : "qqHH_CV_1p5_C2V_1_kl_1_hbbhtautau",
	 "qqH_htt" : "qqH_htt",
	 "singleT" : "singleT",
	 "ttH_hbb" : "ttH_hbb",
}



unc_sources = {
	 "CMS_scale_es_13TeV_2016_DM0" : "CMS_e_FakeTau_dm0_13TeV",
	 "CMS_scale_es_13TeV_2016_DM1" : "CMS_e_FakeTau_dm1_13TeV",
 	 "CMS_scale_mes_13TeV_2016" : "CMS_m_FakeTau_dm0_13TeV",
 	 "CMS_scale_t_13TeV_2016_DM0" : "CMS_scale_t_dm0_13TeV",
 	 "CMS_scale_t_13TeV_2016_DM1" : "CMS_scale_t_dm1_13TeV",
 	 "CMS_scale_t_13TeV_2016_DM10" : "CMS_scale_t_dm10_13TeV",
 	 "CMS_scale_t_13TeV_2016_DM11" : "CMS_scale_t_dm11_13TeV",
	 "CMS_scale_j_Abs" : "CMS_JES_Abs_13TeV",
	 "CMS_scale_j_Abs_2016" : "CMS_JES_Abs_2016_13TeV",
	 "CMS_scale_j_BBEC1" : "CMS_JES_BBEC1_13TeV",
	 "CMS_scale_j_BBEC1_2016" : "CMS_JES_BBEC1_2016_13TeV",
	 "CMS_scale_j_EC2" : "CMS_JES_EC2_13TeV",
	 "CMS_scale_j_EC2_2016" : "CMS_JES_EC2_2016_13TeV",
	 "CMS_scale_j_FlavQCD" : "CMS_JES_FlavQCD_13TeV",
	 "CMS_scale_j_HF" : "CMS_JES_HF_13TeV",
	 "CMS_scale_j_HF_2016" : "CMS_JES_HF_2016_13TeV",
	 "CMS_scale_j_RelBal" : "CMS_JES_RelBal_13TeV",
	 "CMS_scale_j_RelSample_2016" : "CMS_JES_RelSample_2016_13TeV",
}

with open(args.input, 'r') as f:
    yields = json.load(f)

acceptances = {}
for orig_proc, yield_proc in processes.items():
    if yield_proc not in yields:
        raise RuntimeError('Unable to find yields for {} ({}). Available processes: {}' \
                           .format(orig_proc, yield_proc, ', '.join(sorted(yields.keys()))))
    acceptances[orig_proc] = {}
    central_yield = yields[yield_proc]['Central']
    for orig_unc, yield_unc in unc_sources.items():
        if yield_unc not in yields[yield_proc]:
            available_unc = sorted([ unc for unc in yields[yield_proc].keys() if unc != 'Central' ])
            raise RuntimeError('Unable to find yields for unc_source = {} ({}) for process {} ({}).'
                               ' Available unc_sources: {}' \
                               .format(orig_unc, yield_unc, orig_proc, yield_proc, ', '.join(available_unc)))
        acceptances[orig_proc][orig_unc] = {}
        for unc_scale in [ 'Up', 'Down' ]:
            if unc_scale not in yields[yield_proc][yield_unc]:
                print("Warning: yield for {} {} {} is 0".format(orig_proc, orig_unc, unc_scale))
                shifted_yield = 0
            else:
                shifted_yield = yields[yield_proc][yield_unc][unc_scale]
            acc_corr = shifted_yield / central_yield
            acceptances[orig_proc][orig_unc][unc_scale] = acc_corr

with open(args.output, 'w') as f:
    json.dump(acceptances, f, sort_keys=True, indent=4)
