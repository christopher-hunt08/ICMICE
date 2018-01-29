
"""
  Build an event object from the data event
"""

from analysis import tools
from analysis import hit_types

####################################################################################################
def build_event(event_loader) :
  scifi_event = event_loader("scifi")
  tof_event = event_loader("tof")

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

    temp_track = [None for _ in range(16)]

    for tp in track.scifitrackpoints() :
      plane = tools.calculate_plane_id(1, tp.station(), tp.plane())
      hit = hit_types.AnalysisHit(scifi_track_point=tp)

      temp_track[plane] = hit

    if track.tracker() == 0 :
      analysis_event._AnalysisEvent__tracker0_tracks.append(temp_track)
    elif track.tracker() == 1 :
      analysis_event._AnalysisEvent__tracker1_tracks.append(temp_track)
    else : 
      raise "WTF!?"


  return analysis_event




####################################################################################################
class AnalysisEvent(object) :

  def __init__(self, selection_plane=(0, 1, 0), reference_plane=(1, 0)) :
    self.__tof0_spacepoints = []
    self.__tof1_spacepoints = []
    self.__tof2_spacepoints = []

    self.__tracker0_tracks = []
    self.__tracker1_tracks = []

    self.__selection_plane = tools.calculate_plane_id( *selection_plane )
    if self.__selection_plane > 0 :
      self.__selection_tracker = self.__tracker1_tracks
    else:
      self.__selection_tracker = self.__tracker0_tracks
    self.__selection_plane = abs(self.__selection_plane)

    self.__upstream_ref = tools.calculate_plane_id( 1, *reference_plane )
    self.__downstream_ref = tools.calculate_plane_id( 1, *reference_plane )


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

