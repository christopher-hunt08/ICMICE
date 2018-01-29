

import ROOT
import argparse


class Cut_Base(object) :

  def __init__(self) :
    pass


  def is_cut(self, analysis_event) :
    raise NotImplementedError()


