
from .. import LastAnalysis

from _analysis_base import Analysis_Base

import math
import array


class TrackFindingEfficiency(Analysis_Base) :

  def __init__(self) :
    Analysis_Base.__init__(self, "track_finding_efficiency")

    self.__upstream_expected = 0
    self.__upstream_expected_spacepoints = 0
    self.__upstream_found = 0
    self.__upstream_found_spacepoints = 0
    self.__upstream_superfluous = 0
    self.__downstream_expected = 0
    self.__downstream_expected_spacepoints = 0
    self.__downstream_found = 0
    self.__downstream_found_spacepoints = 0
    self.__downstream_superfluous = 0

    self.__upstream_hard_efficiency = 0.0
    self.__upstream_soft_efficiency = 0.0
    self.__downstream_hard_efficiency = 0.0
    self.__downstream_soft_efficiency = 0.0


  def analyse_event(self, analysis_event, weight=1.0) :
    # Require both TOF spacepoints
    if analysis_event.num_tof1_spacepoints() == 1 \
        and analysis_event.num_tof2_spacepoints() == 1 :

      ## Upstream Requirements
      if analysis_event.num_downstream_tracks() == 1 :
        self.__upstream_expected += 1

#        if analysis_event.num_upstream_spacepoints() >= 4 :
#          self.__upstream_expected_spacepoints += 1
        perfect = True
        for i in range(1, 5) :
          if len(analysis_event.upstream_spacepoints()[i]) != 1 :
            perfect = False
            break
        if perfect :
          self.__upstream_expected_spacepoints += 1
          if analysis_event.num_upstream_tracks() == 1 :
            self.__upstream_found_spacepoints += 1

        if analysis_event.num_upstream_tracks() == 1 :
          self.__upstream_found += 1
          self.__upstream_superfluous += analysis_event.num_upstream_tracks()-1
      else :
        self.__upstream_superfluous += analysis_event.num_upstream_tracks()

      ## Downstream Requirements
      if analysis_event.num_upstream_tracks() == 1 :
        self.__downstream_expected += 1

#        if analysis_event.num_downstream_spacepoints() >= 4 :
#          self.__downstream_expected_spacepoints += 1
        perfect = True
        for i in range(1, 5) :
          if len(analysis_event.downstream_spacepoints()[i]) != 1 :
            perfect = False
            break
        if perfect :
          self.__downstream_expected_spacepoints += 1
          if analysis_event.num_downstream_tracks() == 1 :
            self.__downstream_found_spacepoints += 1

        if analysis_event.num_downstream_tracks() == 1 :
          self.__downstream_found += 1
          self.__downstream_superfluous += analysis_event.num_downstream_tracks()-1
      else :
        self.__downstream_superfluous += analysis_event.num_downstream_tracks()


  def _get_plots(self, plot_dict) :
    pass


  def _get_data(self, data_dict) :
    upstream_dict = {}
    upstream_dict['hard_efficiency'] = self.__upstream_hard_efficiency
    upstream_dict['soft_efficiency'] = self.__upstream_soft_efficiency
    upstream_dict['hard_error'] = self.__upstream_hard_error
    upstream_dict['soft_error'] = self.__upstream_soft_error
    upstream_dict['expected'] = self.__upstream_expected
    upstream_dict['hard_expected'] = self.__upstream_expected_spacepoints
    upstream_dict['found'] = self.__upstream_found
    upstream_dict['superfluous'] = self.__upstream_superfluous

    downstream_dict = {}
    downstream_dict['hard_efficiency'] = self.__downstream_hard_efficiency
    downstream_dict['soft_efficiency'] = self.__downstream_soft_efficiency
    downstream_dict['hard_error'] = self.__downstream_hard_error
    downstream_dict['soft_error'] = self.__downstream_soft_error
    downstream_dict['expected'] = self.__downstream_expected
    downstream_dict['hard_expected'] = self.__downstream_expected_spacepoints
    downstream_dict['found'] = self.__downstream_found
    downstream_dict['superfluous'] = self.__downstream_superfluous

    data_dict['upstream'] = upstream_dict
    data_dict['downstream'] = downstream_dict


  def configure_arguments(self, parser) :
    pass


  def parse_arguments(self, namespace) :
    pass


  def conclude(self) :
    self.__upstream_hard_efficiency = float(self.__upstream_found_spacepoints) / float(self.__upstream_expected_spacepoints)
    self.__upstream_soft_efficiency = float(self.__upstream_found) / float(self.__upstream_expected)
    k = self.__upstream_soft_efficiency
    N = float(self.__upstream_expected_spacepoints)
    self.__upstream_hard_error = (1.0/N)*math.sqrt(k*(1.0-k/N))
    N = float(self.__upstream_expected)
    self.__upstream_soft_error = (1.0/N)*math.sqrt(k*(1.0-k/N))

    self.__downstream_hard_efficiency = float(self.__downstream_found_spacepoints) / float(self.__downstream_expected_spacepoints)
    self.__downstream_soft_efficiency = float(self.__downstream_found) / float(self.__downstream_expected)
    k = self.__downstream_soft_efficiency
    N = float(self.__downstream_expected_spacepoints)
    self.__downstream_hard_error = (1.0/N)*math.sqrt(k*(1.0-k/N))
    N = float(self.__downstream_expected)
    self.__downstream_soft_error = (1.0/N)*math.sqrt(k*(1.0-k/N))

