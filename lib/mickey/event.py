
"""
  Build an event object from the data event
"""

from analysis import tools
from analysis import hit_types
from analysis import analysis_track

####################################################################################################
def build_event(event_loader, mc_lookup=None) :
  scifi_event = event_loader("scifi")
  tof_event = event_loader("tof")
  global_event = event_loader("global")

  analysis_event = AnalysisEvent()

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

    trackpoints = [None for _ in range(16)]

    for tp in track.scifitrackpoints() :
      plane = tools.calculate_plane_id(1, tp.station(), tp.plane())
      hit = hit_types.AnalysisHit(scifi_track_point=tp)

      trackpoints[plane] = hit

    temp_track = analysis_track.AnalysisTrack(trackpoints, chisq=track.chi2(), ndf=track.ndf(), p_value=track.P_value())

    if track.tracker() == 0 :
      analysis_event._AnalysisEvent__tracker0_tracks.append(temp_track)
    elif track.tracker() == 1 :
      analysis_event._AnalysisEvent__tracker1_tracks.append(temp_track)
    else : 
      raise "WTF!?"


  global_tracks = global_event.get_tracks()

  for track in global_tracks :
    track.SortTrackPointsByZ()
    trackpoints = []

    for tp in track.GetTrackPoints() :
      hit = hit_types.AnalysisHit(global_track_point=tp)

      trackpoints.append(hit)

    analysis_event._AnalysisEvent__global_tracks.append( analysis_track.AnalysisTrack(trackpoints) )



  if mc_lookup is not None :
    _fill_mc(event_loader("mc"), analysis_event, mc_lookup)

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
class AnalysisEvent(object) :

  def __init__(self, selection_plane=(0, 1, 0), reference_plane=(1, 0)) :
    self.__tof0_spacepoints = []
    self.__tof1_spacepoints = []
    self.__tof2_spacepoints = []

    self.__tracker0_tracks = []
    self.__tracker1_tracks = []

    self.__global_tracks = []

    self.__mc_trackpoints = {}

    self.__selection_plane = tools.calculate_plane_id( *selection_plane )
    if self.__selection_plane > 0 :
      self.__selection_tracker = self.__tracker1_tracks
    else:
      self.__selection_tracker = self.__tracker0_tracks
    self.__selection_plane = abs(self.__selection_plane)

    self.__upstream_ref = tools.calculate_plane_id( 1, *reference_plane )
    self.__downstream_ref = tools.calculate_plane_id( 1, *reference_plane )


    self.__mc_selection_plane = tools.calculate_plane_id( *selection_plane )
    self.__upstream_mc_ref = tools.calculate_plane_id( 0, *reference_plane )
    self.__downstream_mc_ref = tools.calculate_plane_id( 1, *reference_plane )
    for plane in range(-15, 16) :
      self.__mc_trackpoints[plane] = None


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
    return self.__mc_trackpoints[index]


  def mc_selection_trackpoint(self) :
    return self.__mc_trackpoints[self.__mc_selection_plane]


  def mc_upstream_reference_trackpoint(self) :
    return self.__mc_trackpoints[self.__upstream_mc_ref]


  def mc_downstream_reference_trackpoint(self) :
    return self.__mc_trackpoints[self.__downstream_mc_ref]




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

