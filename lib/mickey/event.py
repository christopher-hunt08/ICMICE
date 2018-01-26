
"""
  Build an event object from the data event
"""

####################################################################################################
def build_event(event_loader) :
  mc_event = event_loader("mc")
  scifi_event = event_loader("scifi")
  tof_event = event_loader("tof")


####################################################################################################
class AnalysisEvent(object) :

  def __init__(self) :
    self.

