import itertools
import math
import os
import yaml

from CombineHarvester.CombineTools.ch import CombineHarvester, SystMap
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.TH1.AddDirectory(False)
ROOT.TH2.AddDirectory(False)
ROOT.gROOT.SetMustClean(False)

def listToVector(x, value_type):
  v = ROOT.std.vector(value_type)(len(x))
  for n in range(len(x)):
    v[n] = x[n]
  return v

def rebinAndFill(new_hist, old_hist):
    epsilon = 1e-7

    def check_range(old_axis, new_axis):
        old_min = old_axis.GetBinLowEdge(1)
        old_max = old_axis.GetBinUpEdge(old_axis.GetNbins())
        new_min = new_axis.GetBinLowEdge(1)
        new_max = new_axis.GetBinUpEdge(new_axis.GetNbins())
        return old_min <= new_min and old_max >= new_max

    def get_new_bin(old_axis, new_axis, bin_id_old):
        old_low_edge = round(old_axis.GetBinLowEdge(bin_id_old), 4)
        old_up_edge = round(old_axis.GetBinUpEdge(bin_id_old), 4)
        bin_low_new = new_axis.FindFixBin(old_low_edge)
        bin_up_new = new_axis.FindFixBin(old_up_edge)

        new_up_edge = new_axis.GetBinUpEdge(bin_low_new)
        if not (bin_low_new == bin_up_new or \
                abs(old_up_edge - new_up_edge) <= epsilon * abs(old_up_edge + new_up_edge) * 2):
            old_bins = [ str(old_axis.GetBinLowEdge(n)) for n in range(1, old_axis.GetNbins() + 2)]
            print('old_bins: [{}]'.format(', '.join(old_bins)))
            raise RuntimeError("Uncompatible bin edges")
        return bin_low_new

    def add_bin_content(bin_old, bin_new):
        cnt_old = old_hist.GetBinContent(bin_old)
        err_old = old_hist.GetBinError(bin_old)
        cnt_new = new_hist.GetBinContent(bin_new)
        err_new = new_hist.GetBinError(bin_new)
        cnt_upd = cnt_new + cnt_old
        err_upd = math.sqrt(err_new ** 2 + err_old ** 2)

        new_hist.SetBinContent(bin_new, cnt_upd)
        new_hist.SetBinError(bin_new, err_upd);

    n_dim = old_hist.GetDimension()
    if new_hist.GetDimension() != n_dim:
        raise RuntimeError("Incompatible number of dimensions")
    if n_dim < 1 or n_dim > 2:
        raise RuntimeError("Unsupported number of dimensions")

    if not check_range(old_hist.GetXaxis(), new_hist.GetXaxis()):
        raise RuntimeError("x ranges are not compatible")

    if n_dim > 1 and not check_range(old_hist.GetYaxis(), new_hist.GetYaxis()):
        raise RuntimeError("y ranges are not compatible")

    for x_bin_old in range(old_hist.GetNbinsX() + 2):
        x_bin_new = get_new_bin(old_hist.GetXaxis(), new_hist.GetXaxis(), x_bin_old)
        if n_dim == 1:
            add_bin_content(x_bin_old, x_bin_new)
        else:
            for y_bin_old in range(old_hist.GetNbinsY() + 2):
                y_bin_new = get_new_bin(old_hist.GetYaxis(), new_hist.GetYaxis(), y_bin_old)
                bin_old = old_hist.GetBin(x_bin_old, y_bin_old)
                bin_new = new_hist.GetBin(x_bin_new, y_bin_new)
                add_bin_content(bin_old, bin_new)

def hasNegativeBins(histogram):
  for n in range(1, histogram.GetNbinsX() + 1):
    if histogram.GetBinContent(n) < 0:
      return True
  return False

class Process:
  def __init__(self, name, hist_name, is_signal=False, is_data=False, is_asimov_data=False):
    self.name = name
    self.hist_name = hist_name
    self.is_signal = is_signal
    self.is_data = is_data
    self.is_asimov_data = is_asimov_data
    self.is_background = not (is_signal or is_data)

    if is_data and is_signal:
      raise RuntimeError("Data and signal flags cannot be set simultaneously")
    if is_asimov_data and not is_data:
      raise RuntimeError("Asimov data flag can only be set for data processes")

    if is_signal:
      self.type = "signal"
    elif is_data:
      self.type = "data"
    else:
      self.type = "background"

  def __str__(self):
    return f"Process({self.name}, type={self.type})"

  @staticmethod
  def extractParameters(name_pattern):
    parameters = []
    idx = 0
    while True:
      pos = name_pattern.find("${", idx)
      if pos == -1:
        break
      end = name_pattern.find("}", pos)
      param_name = name_pattern[pos+2:end]
      parameters.append(param_name)
      idx = end
    return parameters

  @staticmethod
  def applyParameters(name_pattern, parameters, param_values):
    if len(parameters) == 1 and type(param_values) != list:
      param_values = [ param_values ]
    if len(param_values) != len(parameters):
      raise RuntimeError("Invalid parameter entry length")
    for n in range(len(parameters)):
      param_name = parameters[n]
      name_pattern = name_pattern.replace("${" + param_name + "}", str(param_values[n]))
    return name_pattern

  @staticmethod
  def fromConfig(entry):
    if type(entry) == str:
      return [ Process(entry, entry) ]
    if type(entry) != dict:
      raise RuntimeError("Invalid entry type")
    base_name = entry["process"]
    base_hist_name = entry.get("hist_name", base_name)
    is_signal = entry.get("is_signal", False)
    is_data = entry.get("is_data", False)
    is_asimov_data = entry.get("is_asimov_data", False)
    if 'param_values' not in entry:
      return [ Process(base_name, base_hist_name, is_signal=is_signal, is_data=is_data, is_asimov_data=is_asimov_data) ]
    parameters = Process.extractParameters(base_name)
    param_values = entry["param_values"]
    if type(param_values) != list or len(param_values) == 0:
      raise RuntimeError("Invalid parameter values")
    processes = []
    for param_entry in param_values:
      name = Process.applyParameters(base_name, parameters, param_entry)
      hist_name = Process.applyParameters(base_hist_name, parameters, param_entry)
      processes.append(Process(name, hist_name, is_signal=is_signal, is_data=is_data, is_asimov_data=is_asimov_data))
    return processes

class HarvesterInterface:
  def __init__(self, cfg_file, input_path, hist_bins=None):
    self.cb = CombineHarvester()

    self.input_path = input_path
    with open(cfg_file, 'r') as f:
      cfg = yaml.safe_load(f)

    self.analysis = cfg["analysis"]
    self.eras = cfg["eras"]
    self.channels = cfg["channels"]
    self.categories = cfg["categories"]

    self.bins = []
    for era, channel, cat in self.ECC():
      bin = self.getBin(era, channel, cat, return_index=False)
      self.bins.append(bin)

    self.processes = {}
    data_process = None
    has_signal = False
    for process in cfg["processes"]:
      new_processes = Process.fromConfig(process)
      for process in new_processes:
        if process.name in self.processes:
          raise RuntimeError(f"Process name {process.name} already exists")
        print(f"Adding {process}")
        self.processes[process.name] = process
        if process.is_data:
          if data_process is not None:
            raise RuntimeError("Multiple data processes defined")
          data_process = process
        if process.is_signal:
          has_signal = True
    if data_process is None:
      raise RuntimeError("No data process defined")
    if not has_signal:
      raise RuntimeError("No signal process defined")

    self.uncertainties = {}
    for unc_entry in cfg["uncertainties"]:
      unc_name = unc_entry["name"]
      if unc_name in self.uncertainties:
        raise RuntimeError(f"Uncertainty {unc_name} already exists")
      self.uncertainties[unc_name] = unc_entry

    self.hist_bins = hist_bins
    if self.hist_bins is None:
      self.hist_bins = cfg.get("hist_bins", None)
    if self.hist_bins is not None:
      self.hist_bins = listToVector(self.hist_bins, 'float')

    self.input_files = {}
    self.shapes = {}

  def getBin(self, era, channel, category, return_name=True, return_index=True):
    name = f'{era}_{self.analysis}_{channel}_{category}'
    if not return_name and not return_index:
      raise RuntimeError("Invalid argument combination")
    if not return_index:
      return name
    index = self.bins.index(name)
    if not return_name:
      return index
    return (index, name)

  def ECC(self):
    return itertools.product(self.eras, self.channels, self.categories)

  def getInputFileName(self, era):
    file_name = f'{era}.root'
    return os.path.join(self.input_path, file_name)

  def getInputFile(self, era):
    if era not in self.input_files:
      file_name = self.getInputFileName(era)
      file = ROOT.TFile.Open(file_name, "READ")
      if file == None:
        raise RuntimeError(f"Cannot open file {file_name}")
      self.input_files[era] = file
    return self.input_files[era]

  def getShape(self, process, era, channel, category, syst_name=None, syst_scale=None):
    key = (process.name, era, channel, category, syst_name, syst_scale)
    if key not in self.shapes:
      if process.is_asimov_data:
        hist = None
        for bkg_proc in self.processes.values():
          if bkg_proc.is_background:
            bkg_hist = self.getShape(bkg_proc, era, channel, category)
            if hist is None:
              hist = bkg_hist.Clone()
            else:
              hist.Add(bkg_hist)
        if hist is None:
          raise RuntimeError("Cannot create asimov data histogram")
      else:
        file = self.getInputFile(era)
        hist_name = f"{channel}/{category}/{process.hist_name}"
        if syst_name is not None:
          hist_name = f"{hist_name}_{syst_name}{syst_scale}"
        hist = file.Get(hist_name)
        if hist == None:
          raise RuntimeError(f"Cannot find histogram {hist_name} in {file.GetName()}")
      if self.hist_bins is not None:
        new_hist = ROOT.TH1F(hist.GetName(), hist.GetTitle(), len(self.hist_bins) - 1, self.hist_bins.data())
        rebinAndFill(new_hist, hist)
        hist = new_hist
      hist.SetDirectory(0)
      if hasNegativeBins(hist):
        raise RuntimeError(f"Negative bins found in histogram {hist_name}")
      self.shapes[key] = hist
    return self.shapes[key]

  def addProcess(self, proc, era, channel, category):
    bin_idx, bin_name = self.getBin(era, channel, category)
    process = self.processes[proc]
    if process.is_data:
      self.cb.AddObservations(["*"], [self.analysis], [era], [channel], [(bin_idx, bin_name)])
    else:
      self.cb.AddProcesses(["*"], [self.analysis], [era], [channel], [proc], [(bin_idx, bin_name)], process.is_signal)

    shape = self.getShape(process, era, channel, category)
    def setShape(p):
      print(f"Setting shape for {p}")
      p.set_shape(shape, True)
    if process.is_data:
      self.cb.cp().process([proc]).channel([channel]).bin([bin_name]).ForEachObs(setShape)
    else:
      self.cb.cp().process([proc]).channel([channel]).bin([bin_name]).ForEachProc(setShape)

  def addSystematicUncertainty(self, syst_name):
    syst_entry = self.uncertainties[syst_name]
    for era, channel, category in self.ECC():
      bin_idx, bin_name = self.getBin(era, channel, category)
      self.cb.cp().channel([channel]).bin([bin_name]) \
          .AddSyst(self.cb, syst_name, syst_entry["type"], SystMap()(syst_entry["value"]))

  def writeDatacards(self, output):
    os.makedirs(output, exist_ok=True)
    background_names = [ proc_name for proc_name, proc in self.processes.items() if proc.is_background ]
    for proc_name, process in self.processes.items():
      if not process.is_signal: continue
      processes = [proc_name] + background_names
      self.cb.cp().process(processes).WriteDatacard(os.path.join(output, f"datacard_{proc_name}.txt"),
                                                    os.path.join(output, f"{proc_name}.root"))

  def createDatacards(self, output, verbose=1):
    try:
      for era, channel, category in self.ECC():
        for process_name in self.processes.keys():
          self.addProcess(process_name, era, channel, category)
      for unc_name in self.uncertainties.keys():
        self.addSystematicUncertainty(unc_name)
      self.cb.SetAutoMCStats(self.cb, 20, True, 1)
      if verbose > 0:
        self.cb.PrintAll()
      self.writeDatacards(output)
    finally:
      for file in self.input_files.values():
        file.Close()