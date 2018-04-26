
"""
  Build an event object from the data event
"""

from analysis import tools
from analysis import hit_types
from analysis import analysis_track


####################################################################################################
CURRENT_MC_LOOKUP = None


####################################################################################################
def build_event(event_loader, mc_lookup=None, selection_plane=-1, reference_plane=1) :
  scifi_event = event_loader("scifi")
  tof_event = event_loader("tof")
  global_event = event_loader("global")

  global CURRENT_MC_LOOKUP
  CURRENT_MC_LOOKUP = mc_lookup

  analysis_event = AnalysisEvent(selection_plane=selection_plane, reference_plane=reference_plane)

  tof_spacepoints = tof_event.GetTOFEventSpacePoint()

  tof0_sp_size = tof_spacepoints.GetTOF0SpacePointArraySize()
  tof1_sp_size = tof_spacepoints.GetTOF1SpacePointArraySize()
  tof2_sp_size = tof_spacepoints.GetTOF2SpacePointArraySize()

  for tof0_i in range(tof0_sp_size) : 
    analysis_event._AnalysisEvent__tof0_spacepoints.append( hit_types.AnalysisSpacePoint( tof=tof_spacepoints.GetTOF0SpacePointArrayElement(tof0_i) ) )

  for tof1_i in range(tof1_sp_size) : 
    analysis_event._AnalysisEvent__tof1_spacepoints.append( hit_types.AnalysisSpacePoint( tof=tof_spacepoints.GetTOF1SpacePointArrayElement(tof1_i) ) )

  for tof2_i in range(tof2_sp_size) : 
    analysis_event._AnalysisEvent__tof2_spacepoints.append( hit_types.AnalysisSpacePoint( tof=tof_spacepoints.GetTOF2SpacePointArrayElement(tof2_i) ) )


  scifi_tracks = scifi_event.scifitracks()

  for track in scifi_tracks :
    if track.GetAlgorithmUsed() != tools.HELICAL_ALGORITHM_ID :
      continue

    trackpoints = {}

    for tp in track.scifitrackpoints() :
      plane = tools.calculate_plane_id(1, tp.station(), tp.plane())
      hit = hit_types.AnalysisHit(scifi_track_point=tp)

      trackpoints[plane] = hit

    temp_track = analysis_track.AnalysisTrack(trackpoints, chisq=track.chi2(), ndf=track.ndf(), p_value=track.P_value(), status=track.GetWasRefit())

    if track.tracker() == 0 :
      analysis_event._AnalysisEvent__tracker0_tracks.append(temp_track)
    elif track.tracker() == 1 :
      analysis_event._AnalysisEvent__tracker1_tracks.append(temp_track)
    else : 
      raise "WTF!?"


  global_chains = global_event.get_primary_chains()
  for chain in global_chains :
    chain_type = chain.get_chain_type()
    global_tracks = chain.GetMatchedTracks()
    trackpoints = []

    for track in global_tracks :
#      track.SortTrackPointsByZ()

      for tp in track.GetTrackPoints() :
        hit = hit_types.AnalysisHit(global_track_point=tp)

        trackpoints.append(hit)

    analysis_event._AnalysisEvent__global_tracks.append( analysis_track.AnalysisTrack(trackpoints, status=int(chain_type)) )



  if CURRENT_MC_LOOKUP is not None :
#    _fill_mc(event_loader("mc"), analysis_event, mc_lookup)
    _fill_virt(event_loader("mc"), analysis_event)

  return analysis_event


####################################################################################################
def _fill_mc(mc_event, data_event, mc_lookup) :
  virtual_hits_count = mc_event.GetVirtualHitsSize()

  for virt_i in range(virtual_hits_count) :
    virt = mc_event.GetAVirtualHit(virt_i)
    station = virt.GetStationId()

    if station in mc_lookup :
      plane_id = mc_lookup[station]
      hit = hit_types.AnalysisHit(virtual_track_point=virt)
      data_event._AnalysisEvent__mc_trackpoints[plane_id] = hit


####################################################################################################
def _fill_virt(mc_event, data_event) :
  virtual_hits_count = mc_event.GetVirtualHitsSize()

  for virt_i in range(virtual_hits_count) :
    virt = mc_event.GetAVirtualHit(virt_i)

    if virt.GetTrackId() == 1 : # The primary particle only
      station = virt.GetStationId()
      hit = hit_types.AnalysisHit(virtual_track_point=virt)
      data_event._AnalysisEvent__mc_trackpoints[station] = hit


####################################################################################################
class AnalysisEvent(object) :

  def __init__(self, selection_plane=-1, reference_plane=1) :
    self.__tof0_spacepoints = []
    self.__tof1_spacepoints = []
    self.__tof2_spacepoints = []

    self.__tracker0_tracks = []
    self.__tracker1_tracks = []

    self.__global_tracks = []

    self.__mc_trackpoints = {}

    self.__selection_plane = selection_plane
    if self.__selection_plane > 0 :
      self.__selection_tracker = self.__tracker1_tracks
    else:
      self.__selection_tracker = self.__tracker0_tracks
    self.__selection_plane = abs(self.__selection_plane)

    self.__upstream_ref = reference_plane
    self.__downstream_ref = reference_plane

    if CURRENT_MC_LOOKUP is not None :
      self.__mc_selection_plane = CURRENT_MC_LOOKUP[selection_plane]
      self.__upstream_mc_ref = CURRENT_MC_LOOKUP[ -1*reference_plane ] # Upstream planes are negative
      self.__downstream_mc_ref = CURRENT_MC_LOOKUP[ reference_plane ]
    else :
      self.__mc_selection_plane = 0
      self.__upstream_mc_ref = 0
      self.__downstream_mc_ref = 0


  def num_tof0_spacepoints(self) :
    return len(self.__tof0_spacepoints)


  def num_tof1_spacepoints(self) :
    return len(self.__tof1_spacepoints)


  def num_tof2_spacepoints(self) :
    return len(self.__tof2_spacepoints)


  def tof01(self, tof0_index=0, tof1_index=0) :
    return self.__tof1_spacepoints[tof1_index].get_time() - self.__tof0_spacepoints[tof0_index].get_time()


  def tof12(self, tof1_index=0, tof2_index=0) :
    return self.__tof2_spacepoints[tof2_index].get_time() - self.__tof1_spacepoints[tof1_index].get_time()


  def tof0_spacepoint(self, index=0) :
    return self.__tof0_spacepoints[index]


  def tof1_spacepoint(self, index=0) :
    return self.__tof1_spacepoints[index]


  def tof2_spacepoint(self, index=0) :
    return self.__tof2_spacepoints[index]


  def num_upstream_tracks(self) :
    return len(self.__tracker0_tracks)


  def num_downstream_tracks(self) :
    return len(self.__tracker1_tracks)


  def upstream_track(self, track_index=0) :
    return self.__tracker0_tracks[track_index]


  def downstream_track(self, track_index=0) :
    return self.__tracker1_tracks[track_index]


  def selection_trackpoint(self, track_index=0) :
    return self.__selection_tracker[track_index][self.__selection_plane]


  def upstream_reference_trackpoint(self, track_index=0) :
    return self.__tracker0_tracks[track_index][self.__upstream_ref]


  def downstream_reference_trackpoint(self, track_index=0) :
    return self.__tracker1_tracks[track_index][self.__downstream_ref]



  def num_global_tracks(self) :
    return len(self.__global_tracks)


  def global_track(self, track_index=0) :
    return self.__global_tracks[track_index]



  def mc_trackpoint(self, index) :
    plane_id = CURRENT_MC_LOOKUP[index]
    return self.__mc_trackpoints[plane_id]


  def mc_virtual_hit(self, plane_id) :
    return self.__mc_trackpoints[plane_id]


  def mc_selection_trackpoint(self) :
    if self.__mc_selection_plane in self.__mc_trackpoints :
      return self.__mc_trackpoints[self.__upstream_mc_ref]
    else :
      return None


  def mc_upstream_reference_trackpoint(self) :
    if self.__upstream_mc_ref in self.__mc_trackpoints :
      return self.__mc_trackpoints[self.__upstream_mc_ref]
    else :
      return None


  def mc_downstream_reference_trackpoint(self) :
    if self.__downstream_mc_ref in self.__mc_trackpoints :
      return self.__mc_trackpoints[self.__downstream_mc_ref]
    else :
      return None




  def print_me(self) :
    if len(self.__tof0_spacepoints) > 0 and len(self.__tof1_spacepoints) > 0 :
      tof01 = self.tof01()
    else :
      tof01 = -1
    if len(self.__tof1_spacepoints) > 0 and len(self.__tof2_spacepoints) > 0 :
      tof12 = self.tof12()
    else :
      tof12 = -1

    if len(self.__tracker0_tracks) > 0 :
      up_p = self.upstream_reference_trackpoint().get_p()
    else :
      up_p = -1
    if len(self.__tracker1_tracks) > 0 :
      down_p = self.downstream_reference_trackpoint().get_p()
    else :
      down_p = -1

    print "MAUS EVENT"

    print "Selection Plane :", self.__selection_plane
    print "Reference Planes:", self.__upstream_ref, self.__downstream_ref
    print "         TOF0                 TOF1         Upstream SciFi      Downstream SciFi        TOF2         EMR"
    print "---------------------------------------------------------------------------------------------------------"
    print "{0:6d} Spacepoints  {1:6d} Spacepoints   {2:6d} Tracks       {3:6d} Tracks   {4:6d} Spacepoints      NI".format(\
                      len(self.__tof0_spacepoints), len(self.__tof1_spacepoints), len(self.__tracker0_tracks), len(self.__tracker1_tracks), len(self.__tof2_spacepoints) )
    print "p:                                            {0:5.2f}              {1:5.2f}".format(up_p, down_p)
    print "TOFs:            {0:4.2f}                                 {1:4.2f}".format(tof01, tof12)
    print

