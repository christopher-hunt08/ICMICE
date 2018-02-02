
import ROOT
from _cuts_base import Cut_Base


####################################################################################################
class Cut_tof01_spacepoints(Cut_Base) :

  def __init__(self) :
    Cut_Base.__init__(self, "tof01_spacepoints")

    self.__histogram = ROOT.TH2F( "cut_tof01_number_spacepoints", ";TOF0;TOF1", 5, -0.5, 4.5, 5, -0.5, 4.5 )


  def _is_cut(self, analysis_event) :
    if (analysis_event.num_tof0_spacepoints() == 1) and (analysis_event.num_tof1_spacepoints() == 1) :
      return False
    else :
      return True


  def fill_histograms(self, analysis_event) :
    self.__histogram.Fill( analysis_event.num_tof0_spacepoints(), analysis_event.num_tof1_spacepoints() )


  def _get_plots(self, plot_dict) :
    plot_dict["tof01_number_spacepoints"] = self.__histogram


  def _get_data(self, data_dict) :
    pass


  def configure_arguments(self, parser) :
    pass

  def parse_arguments(self, namespace) :
    pass


####################################################################################################
class Cut_tof01_time(Cut_Base) :

  def __init__(self) :
    Cut_Base.__init__(self, "tof01_time")

    self.__histogram = ROOT.TH1F( "cut_tof01_time", ";TOF01  [ns];# Events", 500, 0.0, 100.0 )
    self.__cut_window = [0.0, 100.0]


  def _is_cut(self, analysis_event) :
    if (analysis_event.num_tof0_spacepoints() >= 1) and (analysis_event.num_tof1_spacepoints() >= 1) :
      if analysis_event.tof01() > self.__cut_window[0] and analysis_event.tof01() < self.__cut_window[1] :
        return False

      else :
        return True

    else :
      return True


  def fill_histograms(self, analysis_event) :
    if (analysis_event.num_tof0_spacepoints() >= 1) and (analysis_event.num_tof1_spacepoints() >= 1) :
      self.__histogram.Fill( analysis_event.tof01() )
    else :
      self.__histogram.AddBinContent(0)


  def _get_plots(self, plot_dict) :
    plot_dict["tof01_time"] = self.__histogram


  def _get_data(self, data_dict) :
    pass


  def configure_arguments(self, parser) :
    parser.add_argument( "--tof01_cut", nargs=2, default=[0.0, 100.0], type=float, help="Cut Window for TOF01" )

  def parse_arguments(self, namespace) :
    self.__cut_window = namespace.tof01_cut
