import copy
import math
import re
from enum import Enum

from ..common.tools import importROOT
ROOT = importROOT()

from CombineHarvester.CombineTools.ch import SystMap

class UncertaintyScale(Enum):
  Up = 1
  Down = -1

class UncertaintyType(Enum):
  auto = 0
  lnN = 1
  shape = 2

class Uncertainty:
  def appliesTo(self, process, era, channel, category):
    return  Uncertainty.hasMatch(process, self.processes) \
        and Uncertainty.hasMatch(era, self.eras) \
        and Uncertainty.hasMatch(channel, self.channels) \
        and Uncertainty.hasMatch(category, self.categories)

  def canIgnore(self, value_thr):
    if self.type == UncertaintyType.lnN:
      if type(self.value) == float:
        return abs(self.value) < value_thr
      return abs(self.value[UncertaintyScale.Up]) < value_thr and abs(self.value[UncertaintyScale.Down]) < value_thr

  def valueToMap(self):
    if self.type == UncertaintyType.lnN:
      if type(self.value) == float:
        return SystMap()(1 + self.value)
      v_down = 1 + self.value[UncertaintyScale.Down]
      v_up = 1 + self.value[UncertaintyScale.Up]
      print(f'{self.name}: {v_down} {v_up}')
      return SystMap()((v_down, v_up))
    elif self.type == UncertaintyType.shape:
      return SystMap()(1.0)
    else:
      raise RuntimeError(f"value is not defined for {self.name}")

  def resolveType(self, nominal_shape, shape_variations, auto_lnN_threshold, asym_value_threshold):
    if self.type != UncertaintyType.auto:
      return self

    nominal_shape = nominal_shape.Clone()
    for n in range(0, nominal_shape.GetNbinsX() + 2):
      nominal_shape.SetBinError(n, 0)
    all_compatible = True
    for unc_scale in [ UncertaintyScale.Up, UncertaintyScale.Down ]:
      if not Uncertainty.haveCompatibleShapes(nominal_shape, shape_variations[unc_scale], auto_lnN_threshold):
        all_compatible = False
        break
    unc_copy = copy.deepcopy(self)
    if not all_compatible:
      unc_copy.type = UncertaintyType.shape
    else:
      unc_copy.type = UncertaintyType.lnN
      unc_copy.value = {}
      nominal_integral = nominal_shape.Integral()
      for unc_scale in [ UncertaintyScale.Up, UncertaintyScale.Down ]:
        scaled_integral = shape_variations[unc_scale].Integral()
        unc_copy.value[unc_scale] = (scaled_integral - nominal_integral) / nominal_integral
      pass_zero_check, pass_sign_check, pass_magnitude_check = unc_copy._checkValue(raise_error=False)
      if pass_sign_check and pass_magnitude_check:
        delta = abs(unc_copy.value[UncertaintyScale.Up]) - abs(unc_copy.value[UncertaintyScale.Down])
        if abs(delta) <= asym_value_threshold:
          max_value = max(abs(unc_copy.value[UncertaintyScale.Up]), abs(unc_copy.value[UncertaintyScale.Down]))
          unc_copy.value = math.copysign(max_value, unc_copy.value[UncertaintyScale.Up])
      else:
        del unc_copy.value
        unc_copy.type = UncertaintyType.shape
    return unc_copy

  def _checkValue(self, raise_error=True):
    pass_zero_check = True
    pass_sign_check = True
    pass_magnitude_check = True
    if self.type == UncertaintyType.lnN:
      if type(self.value) == float:
        if self.value == 0:
          pass_zero_check = False
        if abs(self.value) >= 1:
          pass_magnitude_check = False
      elif type(self.value) == dict:
        if self.value[UncertaintyScale.Down] == 0 or self.value[UncertaintyScale.Up] == 0:
          pass_zero_check = False
        if self.value[UncertaintyScale.Down] * self.value[UncertaintyScale.Up] > 0:
          pass_sign_check = False
        if abs(self.value[UncertaintyScale.Down]) >= 1 or abs(self.value[UncertaintyScale.Up]) >= 1:
          pass_magnitude_check = False
      else:
        raise RuntimeError("Invalid lnN uncertainty value")
    if not pass_zero_check and raise_error:
      raise RuntimeError(f"Value of lnN uncertainty {self.name} cannot be zero")
    if not pass_sign_check and raise_error:
      raise RuntimeError(f"Up/down values of lnN uncertainty {self.name} must have opposite signs")
    if not pass_magnitude_check and raise_error:
      raise RuntimeError(f"Value of lnN uncertainty {self.name} must be less than 100%")
    return pass_zero_check, pass_sign_check, pass_magnitude_check

  @staticmethod
  def fitFlat(histogram):
    def constantFunction(x,par):
      return par[0]
    fit_func = ROOT.TF1("fit_func", constantFunction, 0, 10, 1)
    fit_func.SetParameter(0, 1.0)
    histogram.Fit(fit_func, "nq")
    chi2 = fit_func.GetChisquare()
    ndf = fit_func.GetNDF()
    p_value = ROOT.TMath.Prob(chi2, ndf)
    fit_param = fit_func.GetParameter(0)
    fit_param_error = fit_func.GetParError(0)
    return chi2, p_value, fit_param, fit_param_error

  @staticmethod
  def haveCompatibleShapes(hist_a, hist_b, p_thr):
    hist_b = hist_b.Clone()
    hist_b.Divide(hist_a)
    _, p_value, _, _ = Uncertainty.fitFlat(hist_b)
    return p_value > p_thr

  @staticmethod
  def hasMatch(value, patterns):
    if len(patterns) == 0:
      return True
    for pattern in patterns:
      if pattern[0] == '^':
        if re.match(pattern, value):
          return True
      elif value == pattern:
        return True
    return False

  @staticmethod
  def fromConfig(entry):
    unc = Uncertainty()
    unc.name = entry["name"]
    unc.type = UncertaintyType[entry["type"]]
    if unc.type == UncertaintyType.lnN:
      value = entry["value"]
      if type(value) in [ float, int ]:
        unc.value = float(value)
      elif type(value) == list and len(value) == 2:
        unc.value = { UncertaintyScale.Down: value[0], UncertaintyScale.Up: value[1] }
      else:
        raise RuntimeError("Invalid lnN uncertainty value")
      unc._checkValue()

    def getPatternList(key):
      if key not in entry:
        return []
      v = entry[key]
      if type(v) == str:
        return [ v ]
      return v

    unc.processes = getPatternList("processes")
    unc.eras = getPatternList("eras")
    unc.channels = getPatternList("channels")
    unc.categories = getPatternList("categories")

    return unc
