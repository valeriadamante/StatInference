class Process:
  def __init__(self, name, hist_name, is_signal=False, is_data=False, is_asimov_data=False, scale=1):
    self.name = name
    self.hist_name = hist_name
    self.is_signal = is_signal
    self.is_data = is_data
    self.is_asimov_data = is_asimov_data
    self.is_background = not (is_signal or is_data)
    self.scale=scale

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
    scale = entry.get("scale", 1)
    if type(scale) == str:
      scale = eval(scale)
    if 'param_values' not in entry:
      return [ Process(base_name, base_hist_name, is_signal=is_signal, is_data=is_data, is_asimov_data=is_asimov_data,
                       scale=scale) ]
    parameters = Process.extractParameters(base_name)
    param_values = entry["param_values"]
    if type(param_values) != list or len(param_values) == 0:
      raise RuntimeError("Invalid parameter values")
    processes = []
    for param_entry in param_values:
      name = Process.applyParameters(base_name, parameters, param_entry)
      hist_name = Process.applyParameters(base_hist_name, parameters, param_entry)
      processes.append(Process(name, hist_name, is_signal=is_signal, is_data=is_data, is_asimov_data=is_asimov_data,
                               scale=scale))
    return processes