

import ROOT
from _cuts_base import Cut_Base

# chisq ndf
####################################################################################################
class Cut_scifi_upstream_chisq_ndf(Cut_Base) :

  def __init__(self) :
    Cut_Base.__init__(self, "upstream_scifi_chisqndf")
    self.__histogram = ROOT.TH1F( "upstream_scifi_chisqndf", ";#chi^{2} / NDF;# Muons", 500, 0.0, 100.0 )
    self.__chisq_ndf_cut = 1e300


  def _is_cut(self, analysis_event) :
    if analysis_event.num_upstream_tracks() == 0 :
      return True

    if analysis_event.upstream_track().get_chisq_ndf() > self.__chisq_ndf_cut :
      return True

    else :
      return False


  def fill_histograms(self, analysis_event) :
    if analysis_event.num_upstream_tracks() > 0 :
      self.__histogram.Fill( analysis_event.upstream_track().get_chisq_ndf() )


  def _get_plots(self, plot_dict) :
    plot_dict["chisq_ndf"] = self.__histogram


  def _get_data(self, data_dict) :
    pass


  def configure_arguments(self, parser) :
    parser.add_argument( "--upstream_chisq_ndf_cut", default=4.0, type=float, help="Cut on upstream chisq per degree of freedom" )


  def parse_arguments(self, namespace) :
    self.__chisq_ndf_cut = namespace.upstream_chisq_ndf_cut


####################################################################################################
class Cut_scifi_upstream_pt(Cut_Base) :

  def __init__(self) :
    Cut_Base.__init__(self, "upstream_scifi_pt")
    self.__histogram = ROOT.TH1F( "upstream_scifi_pt", ";p_{t}   [MeV/c];# Muons", 400, 0.0, 200.0 )
    self.__pt_cut = 1e300


  def _is_cut(self, analysis_event) :
    if analysis_event.num_upstream_tracks() == 0 :
      return True

    if analysis_event.upstream_reference_trackpoint().get_pt() > self.__pt_cut :
      return True

    else :
      return False


  def fill_histograms(self, analysis_event) :
    if analysis_event.num_upstream_tracks() > 0 :
      self.__histogram.Fill( analysis_event.upstream_reference_trackpoint().get_pt() )


  def _get_plots(self, plot_dict) :
    plot_dict["pt"] = self.__histogram


  def _get_data(self, data_dict) :
    pass


  def configure_arguments(self, parser) :
    parser.add_argument( "--upstream_pt_cut", default=180.0, type=float, help="Cut on upstream transverse momentum" )


  def parse_arguments(self, namespace) :
    self.__pt_cut = namespace.upstream_pt_cut


# pz
# p
# pt
# radius

