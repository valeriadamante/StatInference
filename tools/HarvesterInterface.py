import CombineHarvester.CombineTools.ch as ch
import itertools
import os
import yaml

class HarvesterInterface:
    def __init__(self, cfg_file, input_path):
        self.cb = ch.CombineHarvester()

        self.input_path = input_path
        with open(cfg_file, 'r') as f:
            self.cfg = yaml.safe_load(f)

        self.years = self.cfg["years"]
        self.channels = self.cfg["channels"]

        self.prod_to_cat = {
            "qqHH" : [ "class" + c for c in self.cfg["vbf_classes"] ],
            "ggHH" : self.cfg["ggHH_categories"],
        }

        self.cat_to_prod = { }
        self.categories = []
        for prod, cats in self.prod_to_cat.items():
            for cat in cats:
                if cat in self.cat_to_prod:
                    raise RuntimeError("Duplicated category name")
                self.cat_to_prod[cat] = prod
                self.categories.append(cat)

        self.ana = self.cfg["ana_name"]
        self.era = self.cfg["era"]
        self.mass = self.cfg["mass"]
        self.var_name = self.cfg["var_name"]
        self.ggHH_selection = self.cfg["ggHH_selection"]
        self.vbf_selection = self.cfg["vbf_selection"]

        self.ch_categories = []
        for year, channel, cat in self.YCC():
            ch_cat = self.GetCHCat(year, channel, cat, return_index=False)
            self.ch_categories.append(ch_cat)

        self.signals = [ '{}_{}'.format(s, self.ana) for s in self.cfg["signals"] ]

    def GetCHChannel(self, year, channel):
        return '{}_{}'.format(channel, year)

    def GetCHCat(self, year, channel, category, return_name=True, return_index=True):
        name = "{}_{}_{}_{}_{}".format(self.cat_to_prod[category], self.ana, channel, year, category)
        if not return_name and not return_index:
            raise RuntimeError("Invalid argument combination")
        if not return_index:
            return name
        index = self.ch_categories.index(name)
        if not return_name:
            return index
        return (index, name)

    def YCC(self):
        return itertools.product(self.years, self.channels, self.categories)

    def GetInputFileName(self, year, channel, cat):
        is_vbf = cat in self.prod_to_cat["qqHH"]
        sel_name = self.vbf_selection if is_vbf else self.ggHH_selection
        if len(sel_name):
            sel_name = "_" + sel_name
        file_name = '{}_{}_{}{}.root'.format(year, channel, self.var_name, sel_name)
        return os.path.join(self.input_path, file_name)

    def AddObservations(self, year, channel, category):
        ch_cat = self.GetCHCat(year, channel, category)
        ch_channel = self.GetCHChannel(year, channel)
        self.cb.AddObservations(["*"], [self.ana], [self.era], [ch_channel], [ch_cat])
        file_name = self.GetInputFileName(year, channel, category)
        self.cb.cp().process(['data_obs']).channel([ch_channel]).bin([ch_cat[1]]) \
            .ExtractShapes(file_name, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC")

    def AddProcess(self, proc, year, channel, category):
        ch_cat = self.GetCHCat(year, channel, category)
        ch_channel = self.GetCHChannel(year, channel)
        is_signal = proc in self.signals
        self.cb.AddProcesses(["*"], [self.ana], [self.era], [ch_channel], [proc], [ch_cat], is_signal)
        file_name = self.GetInputFileName(year, channel, category)
        self.cb.cp().process([proc]).channel([ch_channel]).bin([ch_cat[1]]) \
            .ExtractShapes(file_name, "$BIN/$PROCESS", "$BIN/$PROCESS_$SYSTEMATIC")
