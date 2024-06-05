import ROOT
import numpy as np

input_name = "ggRadion_HH_bbWW_M_300.root"
input_file = ROOT.TFile.Open(input_name)

channel_names = [key.GetName() for key in input_file.GetListOfKeys()]

for channel_name in channel_names:
    subdir = input_file.Get(channel_name)
    process_names = [key.GetName() for key in subdir.GetListOfKeys()]

    for process_name in process_names:
        print("Looking at Channel ", channel_name, " and process ", process_name)
        process_hist = subdir.Get(process_name)

        x = np.zeros(process_hist.GetNbinsX())
        x_err2 = np.zeros(process_hist.GetNbinsX())
        for binnum in range(process_hist.GetNbinsX()):
            x[binnum] = process_hist.GetBinContent(binnum + 1)
            x_err2[binnum] = process_hist.GetBinError(binnum + 1)**2

        print("x vals = ", x)
        print("err vals = ", x_err2)
