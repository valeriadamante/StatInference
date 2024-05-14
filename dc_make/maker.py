import itertools
import math
import os
import yaml

from CombineHarvester.CombineTools.ch import CombineHarvester

from StatInference.common.tools import hasNegativeBins, listToVector, rebinAndFill, importROOT
from .process import Process
from .uncertainty import Uncertainty, UncertaintyType, UncertaintyScale
from .model import Model
ROOT = importROOT()

class DatacardMaker:
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

    self.model = Model.fromConfig(cfg["model"])

    self.param_bins = {}
    self.processes = {}
    data_process = None
    has_signal = False
    for process in cfg["processes"]:
      new_processes = Process.fromConfig(process, self.model)
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
          param_bin = self.model.paramStr(process.params)
          if param_bin in self.param_bins:
            raise RuntimeError(f"Signal process with parameters {param_bin} already exists")
          self.param_bins[param_bin] = process.params
    if data_process is None:
      raise RuntimeError("No data process defined")
    if not has_signal:
      raise RuntimeError("No signal process defined")

    self.uncertainties = {}
    for unc_entry in cfg["uncertainties"]:
      unc = Uncertainty.fromConfig(unc_entry)
      if unc.name in self.uncertainties:
        raise RuntimeError(f"Uncertainty {unc.name} already exists")
      self.uncertainties[unc.name] = unc

    self.autoMCStats = cfg.get("autoMCStats", { 'apply': False })

    self.hist_bins = hist_bins
    if self.hist_bins is None:
      self.hist_bins = cfg.get("hist_bins", None)
    if self.hist_bins is not None:
      self.hist_bins = listToVector(self.hist_bins, 'double')

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

  def cbCopy(self, param_str, process, era, channel, category):
    bin_idx, bin_name = self.getBin(era, channel, category)
    return self.cb.cp().mass([param_str]).process([process]).bin([bin_name])

  def ECC(self):
    return itertools.product(self.eras, self.channels, self.categories)

  def PPECC(self):
    param_bins = self.param_bins.keys() if self.model.param_dependent_bkg else [ "*" ]
    return itertools.product(self.processes.keys(), param_bins, self.eras, self.channels, self.categories)

  def getInputFile(self, era, model_params):
    file_name = self.model.getInputFileName(era, model_params)
    if file_name not in self.input_files:
      full_file_name = os.path.join(self.input_path, file_name)
      file = ROOT.TFile.Open(full_file_name, "READ")
      if file == None:
        raise RuntimeError(f"Cannot open file {full_file_name}")
      self.input_files[file_name] = file
    return self.input_files[file_name]

  def getShape(self, process, era, channel, category, model_params, unc_name=None, unc_scale=None):
    key = (process.name, era, channel, category, unc_name, unc_scale)
    if key not in self.shapes:
      if process.is_data and (unc_name is not None or unc_scale is not None):
        raise RuntimeError("Cannot apply uncertainty to the data process")
      if process.is_asimov_data:
        hist = None
        for bkg_proc in self.processes.values():
          if bkg_proc.is_background:
            bkg_hist = self.getShape(bkg_proc, era, channel, category, model_params)
            if hist is None:
              hist = bkg_hist.Clone()
            else:
              hist.Add(bkg_hist)
        if hist is None:
          raise RuntimeError("Cannot create asimov data histogram")
      else:
        file = self.getInputFile(era, model_params)
        hist_name = f"{channel}/{category}/{process.hist_name}"
        if unc_name is not None:
          hist_name = f"{hist_name}_{unc_name}{unc_scale}"
        hist = file.Get(hist_name)
        if hist == None:
          raise RuntimeError(f"Cannot find histogram {hist_name} in {file.GetName()}")
      if self.hist_bins is not None:
        new_hist = ROOT.TH1F(hist.GetName(), hist.GetTitle(), len(self.hist_bins) - 1, self.hist_bins.data())
        rebinAndFill(new_hist, hist)
        hist = new_hist
      hist.SetDirectory(0)
      if process.scale != 1:
        hist.Scale(process.scale)
      if hasNegativeBins(hist):
        axis = hist.GetXaxis()
        bins_edges = [ str(axis.GetBinLowEdge(n)) for n in range(1, axis.GetNbins() + 2)]
        bin_values = [ str(hist.GetBinContent(n)) for n in range(1, axis.GetNbins() + 1)]
        print(f'bins_edges: [ {", ".join(bins_edges)} ]')
        print(f'bin_values: [ {", ".join(bin_values)} ]')
        raise RuntimeError(f"Negative bins found in histogram {hist_name}")
      self.shapes[key] = hist
    return self.shapes[key]


  def addProcess(self, proc, era, channel, category):
    bin_idx, bin_name = self.getBin(era, channel, category)
    process = self.processes[proc]

    def add(model_params, param_str):
      if process.is_data:
        self.cb.AddObservations([param_str], [self.analysis], [era], [channel], [(bin_idx, bin_name)])
      else:
        self.cb.AddProcesses([param_str], [self.analysis], [era], [channel], [proc], [(bin_idx, bin_name)], process.is_signal)
      shape = self.getShape(process, era, channel, category, model_params)
      shape_set = False
      def setShape(p):
        nonlocal shape_set
        print(f"Setting shape for {p}")
        if shape_set:
          raise RuntimeError("Shape already set")
        p.set_shape(shape, True)
        shape_set = True
      cb_copy = self.cbCopy(param_str, proc, era, channel, category)
      if process.is_data:
        cb_copy.ForEachObs(setShape)
      else:
        cb_copy.ForEachProc(setShape)


    if process.is_signal:
      model_params = process.params
      param_str = self.model.paramStr(model_params)
      add(model_params, param_str)
    elif self.model.param_dependent_bkg:
      for signal_proc in self.processes.values():
        if signal_proc.is_signal:
          model_params = signal_proc.params
          param_str = self.model.paramStr(model_params)
          add(model_params, param_str)
    else:
      add(None, "*")

  def addUncertainty(self, unc_name):
    unc = self.uncertainties[unc_name]
    for proc, param_str, era, channel, category in self.PPECC():
      process = self.processes[proc]
      if process.is_data: continue
      if not unc.appliesTo(proc, era, channel, category): continue
      model_params = self.param_bins.get(param_str, None)
      if not process.hasCompatibleModelParams(model_params): continue

      systMap = unc.valueToMap()
      cb_copy = self.cbCopy(param_str, proc, era, channel, category)
      cb_copy.AddSyst(self.cb, unc_name, unc.type.name, systMap)
      if unc.type == UncertaintyType.shape:
        shapes = {}
        model_params = self.param_bins.get(param_str, None)
        nominal_shape = self.getShape(self.processes[proc], era, channel, category, model_params)
        for unc_scale in [ UncertaintyScale.Up, UncertaintyScale.Down ]:
          shapes[unc_scale] = self.getShape(self.processes[proc], era, channel, category, model_params,
                                            unc_name, unc_scale.name)
        shape_set = False
        def setShape(syst):
          nonlocal shape_set
          print(f"Setting unc shape for {syst}")
          if shape_set:
            raise RuntimeError("Shape already set")
          syst.set_shapes(shapes[UncertaintyScale.Up], shapes[UncertaintyScale.Down], nominal_shape)
          shape_set = True
        cb_copy = self.cbCopy(param_str, proc, era, channel, category).syst_name([unc_name])
        cb_copy.ForEachSyst(setShape)

  def writeDatacards(self, output):
    os.makedirs(output, exist_ok=True)
    background_names = [ proc_name for proc_name, proc in self.processes.items() if proc.is_background ]
    for proc_name, process in self.processes.items():
      if not process.is_signal: continue
      processes = [proc_name] + background_names
      param_str = self.model.paramStr(process.params)
      dc_file = os.path.join(output, f"datacard_{proc_name}.txt")
      shape_file = os.path.join(output, f"{proc_name}.root")
      self.cb.cp().mass([param_str]).process(processes).WriteDatacard(dc_file, shape_file)

  def createDatacards(self, output, verbose=1):
    try:
      for era, channel, category in self.ECC():
        for process_name in self.processes.keys():
          self.addProcess(process_name, era, channel, category)
      for unc_name in self.uncertainties.keys():
        self.addUncertainty(unc_name)
      if self.autoMCStats["apply"]:
        self.cb.SetAutoMCStats(self.cb, self.autoMCStats["threshold"], self.autoMCStats["apply_to_signal"],
                               self.autoMCStats["mode"])
      if verbose > 0:
        self.cb.PrintAll()
      self.writeDatacards(output)
    finally:
      for file in self.input_files.values():
        file.Close()