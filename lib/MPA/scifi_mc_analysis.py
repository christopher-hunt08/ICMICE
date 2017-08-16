#!/usr/bin/env python

# This file is part of MAUS: http://micewww.pp.rl.ac.uk/projects/maus
#
# MAUS is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MAUS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with MAUS.  If not, see <http://www.gnu.org/licenses/>.
#

# pylint: disable = W0311, E1101, W0613, C0111, R0911, W0621, C0103, R0902

import ROOT
import os
import math
import array
import numpy

import framework
import scifi_extractors
import analysis
from analysis import tools
from analysis import covariances
from analysis import hit_types

"""
  SciFi MC Analysis Classes are stored here.
"""

################################################################################


class scifi_mcanalyser(framework.processor_base) :
  def __init__(self, 'scifi_mc_residuals') :
    pass

  def get_dependencies(inserter) :
    pass

  def _reset(self) :
    pass

  def _process(self, file_reader) :
    pass

  def _store_plots(self, plot_dict) :
    pass

  def _store_data(self, data_dict) :
    pass


