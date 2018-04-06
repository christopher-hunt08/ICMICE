

import ROOT
from _cuts_base import Cut_Base

# chisq ndf
####################################################################################################
class SciFiUpstreamChisqNDF(Cut_Base) :

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
    parser.add_argument( "--upstream_chisq_ndf_cut", default=100.0, type=float, help="Cut on upstream chisq per degree of freedom" )


  def parse_arguments(self, namespace) :
    self.__chisq_ndf_cut = namespace.upstream_chisq_ndf_cut


####################################################################################################
class SciFiUpstreamPt(Cut_Base) :

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
    parser.add_argument( "--upstream_pt_cut", default=10000.0, type=float, help="Cut on upstream transverse momentum" )


  def parse_arguments(self, namespace) :
    self.__pt_cut = namespace.upstream_pt_cut


####################################################################################################
class SciFiRefitStatus(Cut_Base) :

  def __init__(self) :
    Cut_Base.__init__(self, "scifi_refit_status")
    self.__cut_refits = False
    self.__histogram_upstream = ROOT.TH1F( "upstream_refit_status", ";Status;# Muons", 5, -0.5, 4.5 )
    self.__histogram_downstream = ROOT.TH1F( "downstream_refit_status", ";Status;# Muons", 5, -0.5, 4.5 )


  def _is_cut(self, analysis_event) :
    if not self.__cut_refits :
      return False

    if analysis_event.num_upstream_tracks() > 0 :
      if analysis_event.upstream_track().get_status() > 0 :
        return True

    if analysis_event.num_downstream_tracks() > 0 :
      if analysis_event.downstream_track().get_status() > 0 :
        return True

    else :
      return False


  def fill_histograms(self, analysis_event) :
    if analysis_event.num_upstream_tracks() > 0 :
      self.__histogram_upstream.Fill( analysis_event.upstream_track().get_status() )
    if analysis_event.num_downstream_tracks() > 0 :
      self.__histogram_downstream.Fill( analysis_event.downstream_track().get_status() )


  def _get_plots(self, plot_dict) :
    plot_dict["refit_status_upstream"] = self.__histogram_upstream
    plot_dict["refit_status_downstream"] = self.__histogram_downstream


  def _get_data(self, data_dict) :
    pass


  def configure_arguments(self, parser) :
    parser.add_argument( "--scifi_refit_cut", action="store_true", help="Throw tracks that have been refitted by some means." )


  def parse_arguments(self, namespace) :
    self.__cut_refits = namespace.scifi_refit_cut


####################################################################################################
class SciFiTransmission(Cut_Base) :

  def __init__(self) :
    Cut_Base.__init__(self, "scifi_transmission_cut")
    self.__histogram = ROOT.TH2F( "scifi_track_numbers", ";# Upstream Tracks;# Downstream Tracks", 5, -0.5, 4.5, 5, -0.5, 4.5 )
    self.__do_cut = False


  def _is_cut(self, analysis_event) :
    if not self.__do_cut :
      return False

    ups = analysis_event.num_upstream_tracks()
    downs = analysis_event.num_downstream_tracks()

    if ups == 1 and downs == 1 : 
      return False
    else :
      return True



  def fill_histograms(self, analysis_event) :
    ups = analysis_event.num_upstream_tracks()
    downs = analysis_event.num_downstream_tracks()

    self.__histogram.Fill(ups, downs)


  def _get_plots(self, plot_dict) :
    plot_dict["scifi_track_numbers"] = self.__histogram


  def _get_data(self, data_dict) :
    pass


  def configure_arguments(self, parser) :
    parser.add_argument( "--scifi_transmission", action="store_true", help="Only analyse events that have precisely 1 track in each tracker" )


  def parse_arguments(self, namespace) :
    self.__do_cut = namespace.scifi_transmission





