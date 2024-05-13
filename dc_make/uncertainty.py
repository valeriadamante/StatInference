from enum import Enum
import re

from CombineHarvester.CombineTools.ch import SystMap

class UncertaintyScale(Enum):
  Up = 1
  Down = -1

class UncertaintyType(Enum):
  auto = 0
  lnN = 1
  shape = 2

class Uncertainty:
  def __init__(self):
    pass

  def appliesTo(self, process, era, channel, category):
    return  Uncertainty.hasMatch(process, self.processes) \
        and Uncertainty.hasMatch(era, self.eras) \
        and Uncertainty.hasMatch(channel, self.channels) \
        and Uncertainty.hasMatch(category, self.categories)

  def valueToMap(self):
    if self.type == UncertaintyType.lnN:
      if type(self.value) == float:
        return SystMap()(1 + self.value)
      v_down = 1 - self.value[UncertaintyScale.Down]
      v_up = 1 + self.value[UncertaintyScale.Up]
      return SystMap()(v_down, v_up)
    elif self.type == UncertaintyType.shape:
      return SystMap()(1.0)
    else:
      raise RuntimeError(f"value is not defined for {self.name}")

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
