import sys
import ROOT

file_name = sys.argv[1]

ws_file = ROOT.TFile(file_name, 'READ')
ws   = ws_file.Get("w")
data = ws.data("data_obs")
mc_model = ws.genobj("ModelConfig")

nll_args = ROOT.RooLinkedList()
nll_args.Add(ROOT.RooFit.Constrain(mc_model.GetNuisanceParameters()))
nll_args.Add(ROOT.RooFit.Extended(mc_model.GetPdf().canBeExtended()))
nll = mc_model.GetPdf().createNLL(data, nll_args)

if len(sys.argv) > 2:
    snapshot = sys.argv[2]
    ws.loadSnapshot(snapshot)

print(nll.getVal())
