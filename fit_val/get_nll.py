import os
import sys
import ROOT

file_name = sys.argv[1]

print("Loading workspace...")
ws_file = ROOT.TFile(file_name, 'READ')
ws   = ws_file.Get("w")
data = ws.data("data_obs")
mc_model = ws.genobj("ModelConfig")

print("Creating NLL...")
nll_args = ROOT.RooLinkedList()
nll_args.Add(ROOT.RooFit.Constrain(mc_model.GetNuisanceParameters()))
nll_args.Add(ROOT.RooFit.Extended(mc_model.GetPdf().canBeExtended()))
nll = mc_model.GetPdf().createNLL(data, nll_args)

if len(sys.argv) > 2:
    p = sys.argv[2]
    if os.path.isfile(p):
        print("Loading fit results...")
        fit_file = ROOT.TFile(p, 'READ')
        bestfit = fit_file.Get('fit_s')
        fit_args= bestfit.floatParsFinal()
        snapshot_name = "bestfitparams"
        ws.saveSnapshot(snapshot_name, ROOT.RooArgSet(fit_args), True)
    else:
        snapshot_name = p
    print("Loading snapshot...")
    ws.loadSnapshot(snapshot_name)

print('NLL = {}'.format(nll.getVal()))
