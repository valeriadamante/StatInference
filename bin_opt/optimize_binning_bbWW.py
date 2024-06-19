import ROOT
import numpy as np
import os

input_name = "ggRadion_HH_bbWW_M_300.root"
#input_name = "/eos/user/d/daebi/shape_files_15Apr2022/shape_m300.root"
#input_name = "shape_m300.root"
input_file = ROOT.TFile.Open(input_name)

channel_names = [key.GetName() for key in input_file.GetListOfKeys()]



for channel_name in channel_names:
    print("At channel ", channel_name)

    hist_dict = {}
    hist_dict['total'] = ROOT.TH1D('total', 'total', 1000, 0.0, 1.0)
    subdir = input_file.Get(channel_name)
    process_names = [key.GetName() for key in subdir.GetListOfKeys()]

    for process_name in process_names:
        if process_name == 'data_obs': continue
        #print("Looking at Channel ", channel_name, " and process ", process_name)
        process_hist = subdir.Get(process_name)

        x = np.zeros(process_hist.GetNbinsX())
        x_err2 = np.zeros(process_hist.GetNbinsX())
        for binnum in range(process_hist.GetNbinsX()):
            x[binnum] = process_hist.GetBinContent(binnum + 1)
            x_err2[binnum] = process_hist.GetBinError(binnum + 1)**2

        #print("x vals = ", x)
        #print("err vals = ", x_err2)

        hist_dict[process_name] = process_hist
        #print("Process ", process_name, ' has bins ', process_hist.GetNbinsX())
        hist_dict['total'].Add(process_hist)



    tmp_signal = hist_dict['ggRadion_HH_bbWW_M_300'].Integral()
    nbins = 4

    bin_content_goal = tmp_signal/nbins

    print("Channel ", channel_name, " has a total of ", tmp_signal, " entries, so the per-bin goal is ", bin_content_goal)

    custom_binning = [1.0]
    custom_totals = []
    done = False
    start = hist_dict['ggRadion_HH_bbWW_M_300'].GetNbinsX()+1
    while not done:
        #for bincount in range(start, hist_dict['ggRadion_HH_bbWW_M_300'].GetNbinsX()+1):
        for bincount in range(start, -1, -1):
            if hist_dict['ggRadion_HH_bbWW_M_300'].Integral(bincount, start) >= bin_content_goal:
                custom_binning.append(hist_dict['ggRadion_HH_bbWW_M_300'].GetXaxis().GetBinLowEdge(bincount))
                custom_totals.append(hist_dict['ggRadion_HH_bbWW_M_300'].Integral(bincount, start))
                print("At bin ", bincount, " integral is ", hist_dict['ggRadion_HH_bbWW_M_300'].Integral(bincount, start))
                start = bincount-1
                break
            if bincount ==  0:
                custom_binning.append(bincount)
                custom_totals.append(hist_dict['ggRadion_HH_bbWW_M_300'].Integral(bincount, start))
                done = True
    print("Found the bin edge set, it is ", custom_binning)
    print("With totals per bin as ", custom_totals)

    custom_binning.sort()


    #Now create a new yaml for just this channel and save it somewhere
    yaml_example = "../config/x_hh_bbww_run3.yaml"
    if not os.path.exists("tmp_yamls/"):
        os.makedirs("tmp_yamls/")
    yaml_name = channel_name
    new_yaml = open("tmp_yamls/"+yaml_name+".yaml", 'w')
    old_yaml = open(yaml_example, 'r')

    for line in old_yaml:
        if 'ee' in line:
            if 'ee' in channel_name.split('_'):
                new_yaml.write(line)
        elif 'emu' in line:
            if 'emu' in channel_name.split('_'):
                new_yaml.write(line)
        elif 'mumu' in line:
            if 'mumu' in channel_name.split('_'):
                new_yaml.write(line)
        elif 'res1b' in line:
            if 'res1b' in channel_name.split('_'):
                new_yaml.write(line)
        elif 'res2b' in line:
            if 'res2b' in channel_name.split('_'):
                new_yaml.write(line)
        elif 'boost' in line:
            if 'boost' in channel_name.split('_'):
                new_yaml.write(line)
        elif 'hist_bins:' in line:
            new_yaml.write("hist_bins: {custom_binning}\n".format(custom_binning = custom_binning))
        else:
            new_yaml.write(line)
    print("Created new yaml"+yaml_name)

    #if not os.path.exists("tmp_datacards/"):
    #    os.makedirs("tmp_datacards/")

    #os.system("cmbEnv python3 StatInference/dc_make/create_datacards.py --input /eos/user/d/daebi/bbWW_shape_files --output StatInference/bin_opt/tmp_shapes --config StatInference/bin_opt/tmp_datacards/Run3_2022_hh_bbww_mumu_res2b.yaml")