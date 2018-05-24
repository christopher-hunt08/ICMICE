
from .. import LastAnalysis

from _analysis_base import Analysis_Base

import ROOT
import math
import array


class TOFAnalysis(Analysis_Base) :

  def __init__(self) :
    Analysis_Base.__init__(self, "TOF_analysis")

    self.__histogram = ROOT.TH2F("tof01_12", ";TOF01;TOF12", 500, 0.0, 100.0, 500, 0.0, 100.0)


  def analyse_event(self, analysis_event, weight=1.0) :
    # Require all TOF spacepoints
    if analysis_event.num_tof0_spacepoints() != 1 :
      return
    if analysis_event.num_tof1_spacepoints() != 1 :
      return
    if analysis_event.num_tof2_spacepoints() != 1 :
      return

    self.__histogram.Fill(analysis_event.tof01(), analysis_event.tof12())


  def _get_plots(self, plot_dict) :
    plot_dict['tof01-12'] = self.__histogram


  def _get_data(self, data_dict) :
    pass


  def configure_arguments(self, parser) :
    pass


  def parse_arguments(self, namespace) :
    pass


  def conclude(self) :
    pass

