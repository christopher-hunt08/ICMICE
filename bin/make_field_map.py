#!/usr/bin/env python

#  This file is part of MAUS: http://micewww.pp.rl.ac.uk/projects/maus
#
#  MAUS is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  MAUS is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with MAUS.  If not, see <http://www.gnu.org/licenses/>.

"""
Example to demonstrate how to make a field map.
"""

import Configuration  # MAUS configuration (datacards) module
import maus_cpp.globals as maus_globals # MAUS C++ calls
import maus_cpp.field as field # MAUS field map calls
import xboa.Common  # xboa post-processor library
import math
import ROOT

ROOT.gROOT.SetBatch(True)

RESOLUTION = 50.0
ZERO_FIELD = 1.0e-4

def main():
    """
    Make a plot of z, bz
    """
    # set up datacards
    print "Welcome to MAUS field map maker example"
    configuration = Configuration.Configuration().\
                                          getConfigJSON(command_line_args=True)

    # initialise field maps and geometry
    print "Building field maps and other initialisation (this can take a while)"
    maus_globals.birth(configuration)
    geometry = maus_globals.get_monte_carlo_mice_modules()

    # Find the size of the geomtry
    length = geometry.get_property("Dimensions", "hep3vector")["z"]

    # make lists of z, bz points
    print "Getting field values"

    z_list = []
    bz_list = []

    z = -0.5*length
    max_z = 0.5*length
    trim_low = 0
    trim_high = 0
    for i in range( int(length / RESOLUTION) ) :
        if z > max_z :
            z = max_z

        (bx_field, by_field, bz_field, ex_field, ey_field, ez_field) = \
                                        field.get_field_value(0., 0., z, 0.)

        B = math.sqrt( bx_field**2 + by_field**2 + bz_field**2 )

        if B < ZERO_FIELD :
            if trim_high == trim_low :
                trim_low = i
                trim_high = i
        else :
            trim_high = i

        z_list.append(z)  # z in millimetres
        bz_list.append(bz_field*1e3)  # bz in T
        print 'z:', z, ' ** b:', bx_field, by_field, bz_field, \
                               'e:', ex_field, ey_field, ez_field
        z += RESOLUTION

    # Cut off the areas with no field
    print "Trimming Tails"

    z_list = z_list[trim_low : trim_high]
    bz_list = bz_list[trim_low : trim_high]

    # now make a ROOT graph of bz against z
    print "Graphing field values"
    canvas = xboa.Common.make_root_canvas("bz vs z")
    hist, graph = xboa.Common.make_root_graph("bz vs z", z_list, "z [m]",
                                                         bz_list, "B_{z} [T]")
    hist.Draw()
    graph.Draw('l')
    canvas.Update()
    canvas.Print('bfield_vs_z.pdf')

    outfile = ROOT.TFile("bfield_vs_z.root", "RECREATE")
    graph.Write("bfield_vs_z")
    outfile.Close()

    # Clean up
    maus_globals.death()
    # Finished
    print "Finished"

if __name__ == "__main__":
    main()
