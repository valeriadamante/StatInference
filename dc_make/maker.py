import itertools
import math
import os
import yaml

from CombineHarvester.CombineTools.ch import CombineHarvester

from StatInference.common.tools import hasNegativeBins, listToVector, rebinAndFill, importROOT
from .process import Process
from .uncertainty import Uncertainty, UncertaintyType, UncertaintyScale
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
      unc = Uncertainty.fromConfig(unc_entry)
      if unc.name in self.uncertainties:
        raise RuntimeError(f"Uncertainty {unc.name} already exists")
      self.uncertainties[unc.name] = unc

    self.autoMCStats = cfg.get("autoMCStats", { 'apply': False })

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

  def cbCopy(self, process, era, channel, category):
    bin_idx, bin_name = self.getBin(era, channel, category)
    return self.cb.cp().process([process]).bin([bin_name])

  def ECC(self):
    return itertools.product(self.eras, self.channels, self.categories)

  def PECC(self):
    return itertools.product(self.processes.keys(), self.eras, self.channels, self.categories)

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
      if process.scale != 1:
        hist.Scale(process.scale)
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
    cb_copy = self.cbCopy(proc, era, channel, category)
    if process.is_data:
      cb_copy.ForEachObs(setShape)
    else:
      cb_copy.ForEachProc(setShape)

  def addUncertainty(self, unc_name):
    unc = self.uncertainties[unc_name]
    for proc, era, channel, category in self.PECC():
      if unc.appliesTo(proc, era, channel, category):
        systMap = unc.valueToMap()
        cb_copy = self.cbCopy(proc, era, channel, category)
        cb_copy.AddSyst(self.cb, unc_name, unc.type.name, systMap)
        if unc.type == UncertaintyType.shape:
          shapes = {}
          nominal_shape = self.getShape(self.processes[proc], era, channel, category)
          for unc_scale in [ UncertaintyScale.Up, UncertaintyScale.Down ]:
            shapes[unc_scale] = self.getShape(self.processes[proc], era, channel, category, unc_name, unc_scale.name)
          def setShape(syst):
            print(f"Setting unc shape for {syst}")
            syst.set_shapes(shapes[UncertaintyScale.Up], shapes[UncertaintyScale.Down], nominal_shape)
          cb_copy = self.cbCopy(proc, era, channel, category).syst_name([unc_name])
          cb_copy.ForEachSyst(setShape)

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