import argparse
import re

parser = argparse.ArgumentParser(description='Plot shape variations.')
parser.add_argument('--input', required=True, type=str, help="input root file")
parser.add_argument('--output', required=True, type=str, help="output pdf file")
parser.add_argument('--bin', required=False, default=None, type=str, help="bin")
parser.add_argument('--process', required=True, type=str, help="process")
parser.add_argument('--unc', required=True, type=str, help="uncertainty")
parser.add_argument('--legend-pos', required=False, default='right', type=str, help="uncertainty")
args = parser.parse_args()

import ROOT
ROOT.gROOT.SetBatch(True)

def ToFixedBinSize(hist):
    hist_fixed = ROOT.TH1F('', '', hist.GetNbinsX(), 0.5, hist.GetNbinsX() + 0.5)
    hist_fixed.SetDirectory(0)
    for n in range(0, hist.GetNbinsX() + 2):
        hist_fixed.SetBinContent(n, hist.GetBinContent(n))
        hist_fixed.SetBinError(n, hist.GetBinError(n))
    return hist_fixed

def GetRelShift(central, shifted):
    rel = shifted.Clone()
    rel.Add(central, -1)
    rel.Divide(central)
    return rel, (shifted.Integral() - central.Integral()) / central.Integral()

def GetShapes(input_file, bin, process_name, unc):
    hists = {}
    shape_names = [ 'central', 'Up', 'Down' ]
    input_root = ROOT.TFile.Open(input_file)
    if input_root == None:
        raise RuntimeError('{} not found'.format(input_file))
    if bin is None:
        bin_dir = input_root
    else:
        bin_dir = input_root.Get(bin)
        if bin_dir == None:
            available_bins = [ str(key.GetName()) for key in input_root.GetListOfKeys() ]
            raise RuntimeError('{} not found. Available bins: {}'.format(bin, ', '.join(available_bins)))
    hist_names = [ str(key.GetName()) for key in bin_dir.GetListOfKeys() ]

    if process_name == 'bkg':
        #process_regex = re.compile('(?!(data_obs|ggHH.*|qqHH.*|ggGraviton.*|ggRadion.*|DY_[0-9]b.*$))(.*)')
        process_regex = re.compile('(?!(data_obs|ggHH.*|qqHH.*|ggGraviton.*|ggRadion.*|DY$))(.*)')
    else:
        process_regex = re.compile('({})'.format(args.process))
    unc_regex = re.compile('.*(Up|Down)$')

    processes = []
    all_unc = []
    for hist_name in hist_names:
        proc_match = process_regex.match(hist_name)
        unc_match = unc_regex.match(hist_name)
        if proc_match is not None and unc_match is None:
            proc_name = proc_match.groups()[-1]
            if proc_name in processes:
                raise RuntimeError("Duplicated process {}".format(proc_name))
            processes.append(proc_name)
    if len(processes) == 0:
        raise RuntimeError('Process {} not found.'.format(process_name))
    print("Processes: {}".format(', '.join(processes)))
    for proc_name in processes:
        proc_hist_names = [ proc_name, '{}_{}Up'.format(proc_name, unc), '{}_{}Down'.format(proc_name, unc) ]
        for n in range(len(shape_names)):
            hist_name = proc_hist_names[n]
            shape_name = shape_names[n]
            if hist_name not in hist_names:
                print("Warning: {} not found".format(hist_name))
            else:

                hist = bin_dir.Get(hist_name)
                if hist == None:
                    raise RuntimeError("Error while loading {}".format(hist_name))
                hist = hist.Clone()
                hist.SetDirectory(0)
                hist = ToFixedBinSize(hist)
                print("Adding {} to {}. Integral = {}".format(hist_name, shape_name, hist.Integral()))
                if shape_name in hists:
                    hists[shape_name].Add(hist)
                else:
                    hists[shape_name] = hist
        for hist_name in hist_names:
            unc_match = re.match('{}_(.*)(Up|Down)'.format(proc_name), hist_name)
            if unc_match:
                unc_name = unc_match.groups()[0]
                if unc_name not in all_unc:
                    all_unc.append(unc_name)

    input_root.Close()
    shapes = []
    for name in shape_names:
        if name not in hists:
            raise RuntimeError("{} histogram for process {} in bin {} not found. Availalbe uncertainties: {}" \
                               .format(name, process_name, bin, ', '.join(all_unc)))
        print("Total {} integral = {}".format(name, hists[name].Integral()))
        shapes.append(hists[name])
    return tuple(shapes)

central, up, down = GetShapes(args.input, args.bin, args.process, args.unc)

ROOT.gStyle.SetOptStat(0)
c = ROOT.TCanvas('', '', 400, 400)
#c.Draw()
c.cd()
ratio_size = 0.25
main_pad = ROOT.TPad('', '', 0.0, ratio_size, 1., 1.)
ratio_pad = ROOT.TPad('', '', 0.0, 0.0, 1., ratio_size )
left_margin, right_margin = 0.1, 0.05
main_pad.SetLeftMargin(left_margin)
ratio_pad.SetLeftMargin(left_margin)
main_pad.SetRightMargin(right_margin)
ratio_pad.SetRightMargin(right_margin)

main_pad.SetBottomMargin(0)
ratio_pad.SetTopMargin(0)
ratio_pad.SetBottomMargin(0.35)


main_pad.Draw()
ratio_pad.Draw()

axis_title_font_size = 0.05
pad_factor = (1 - ratio_size) / ratio_size

main_pad.cd()
main_pad.SetLogy(1)
bin_name = args.bin if args.bin is not None else ''
central.SetTitle('{} {} {}'.format(bin_name, args.process, args.unc))
central.GetYaxis().SetTitle('Events')
central.GetXaxis().SetTitle('')
central.GetXaxis().SetLabelOffset(10)
central.GetYaxis().SetTitleSize(axis_title_font_size)
central.GetYaxis().SetTitleOffset(1)
central.SetLineColor(ROOT.kBlack)
up.SetLineColor(ROOT.kBlue)
down.SetLineColor(ROOT.kRed)

central.Draw()
up.Draw('SAME')
down.Draw('SAME')

up_central, up_diff = GetRelShift(central, up)
down_central, down_diff = GetRelShift(central, down)

legend_opt = 'le'
if args.legend_pos == 'right':
    legend = ROOT.TLegend(0.6, 0.65, 0.95, 0.9)
else:
    legend = ROOT.TLegend(0.1, 0.65, 0.45, 0.9)
legend.AddEntry(central, 'Nominal ({:.1f})'.format(central.Integral()), legend_opt)
legend.AddEntry(up, 'Up ({:.2f}%)'.format(up_diff * 100), legend_opt)
legend.AddEntry(down, 'Down ({:.2f}%)'.format(down_diff * 100), legend_opt)
legend.Draw()

ratio_pad.cd()
ratio_pad.SetGridy()


up_central.GetYaxis().SetTitle('rel. shift')
up_central.GetYaxis().SetTitleSize(0.1)
up_central.GetYaxis().SetTitleOffset(0.5)
up_central.GetYaxis().SetLabelSize(0.1)
up_central.GetXaxis().SetTitle('dnn bin')
up_central.GetXaxis().SetTitleSize(axis_title_font_size * pad_factor)
up_central.GetXaxis().SetLabelSize(0.13)
up_central.GetXaxis().SetTitleOffset(0)
up_central.SetTitle('')
up_central.Draw()
down_central.Draw('SAME')

main_pad.Modified()
# ratio_pad.Modified()
c.Update()
c.Print(args.output)
