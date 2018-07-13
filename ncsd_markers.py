#!/usr/bin/env python

"""
NCSD markers

a module for Neurolucida Companion for Spatial Distribution Analysis

Analyze the apatial distribution of markers from a XML file produced by Neurolucida.
The following parameters are reported:

- number of markers by contours' sectors
- area of contours' sectors

Copyright Olivier Friard 2010-2011

This file is part of "Neurolucida Companion for Spatial Distribution Analysis".

"Neurolucida Companion for Spatial Distribution Analysis" is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or any later version.

"Neurolucida Companion for Spatial Distribution Analysis" is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with "Neurolucida Companion for Spatial Distribution Analysis"; see the file COPYING.TXT.  If not, write to
the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
"""

from __future__ import division

version = '12'
"""
Module history:
 
- 2011-05-10: improved error handle
- 2011-03-08: added file for EPL count (_epl)
- 2011-02-23: added EPL markers count
              modified file output format
"""

import sys
import os
from coordinate_geometry import *
import xml.etree.ElementTree as ET
from math import *
from svg import *

default_sector_degrees = 20

sy = -1
sector_radius = 2500
colors = {'gran': (255, 228, 0), 'glom': (0, 255, 255), 'epl': (255, 0, 255)}

def check_options():
    """check command line options"""

    from optparse import OptionParser
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-a", "--angle", dest="angle", help="angle of sectors in degree (must be divisor of 360).\nThe default value is 20 degrees")
    parser.add_option("-i", "--input", dest="input", help="Neurolucida XML file(s) (use \"*.xml\" to select more files)")
    parser.add_option("-o", "--output", dest="output", help="write results to file")
    parser.add_option("-v", "--verbose", action="store_true", default=False, dest="verbose", help="Print program messages")
    parser.add_option("-s", "--svg", action="store_true", default=False, dest="svg", help="Write SVG file")

    parser.add_option("-d", "--debug", action="store_true", default=False, dest="debug", help="Add information useful for debugging in output")

    (options, args) = parser.parse_args()

    if options.angle:
        sector_degrees = int(options.angle)
    else:
        sector_degrees = int(default_sector_degrees)

    verbose = options.verbose
    svg = options.svg
    DEBUG = options.debug

    if options.input:
        infile = options.input
        if '*' in infile:
            import glob
            infile = glob.glob(infile)
    else:
        infile = input('Input XML file: ')
        if '*' in infile:
            import glob
            infile = glob.glob(infile)

    if type(infile) == type(''):
        import os
        if not os.path.exists(infile):
            print('\n"%s" not found!\n' % infile)
            sys.exit(1)

    if options.output:
        outfile = options.output
    else:
        outfile = ''
    return (sector_degrees, verbose, svg, infile, outfile, DEBUG)



def main(sector_degrees, verbose, svg, infile, outfile):

    def parse_neurolucida_xml(infile):
        contours = {'glomint':{}, 'glomext':{}, 'gran':{}, 'line':{}}
        markers = {}

        tree = ET.parse(infile)
        root = tree.getroot()

        for contour in root.getiterator('contour'):
            name = contour.findtext('name')
            if verbose:
                print('found "%s" contour' % name)
            for point in contour.getiterator('point'):
                if name in contours:
                    contours[name][len(contours[name])] = (float(point.get('x')), float(point.get('y')))

        for c in ['glomext', 'glomint', 'gran', 'line']:
            if verbose:
                print('"%s" contour vertices number: %d' % (c, len(contours[c])))
            if not len(contours[c]):
                print('\n"%s" contour is invalid (%d vertice(s))\n' % (c, len(contours[c])))
                return None, None

        for marker in root.getiterator('marker'):
            for point in marker.getiterator('point'):
                markers[len(markers)] = (float(point.get('x')), float(point.get('y')))
        if verbose:
            print('Markers number: ', len(markers))

        return contours, markers

    def draw_svg(d, scene, color):

        ### mark 1st point of contour
        #print DEBUG
        #raw_input('key')
        if DEBUG:
            scene.add(Circle((d[0][0], sy * d[0][1]), 5, color, 1))
            # point number

        for p in d:
            n = p + 1
            if n == len(d):
               n = 0

            ### draw point of contour
            if DEBUG:
               scene.add(Circle((d[p][0], sy * d[p][1]), 2, (0, 0, 0), 1))
               ### point number
               scene.add(Text((d[p][0] + 4, sy * d[p][1]), "%d" % (p + 1), 10))

            scene.add(Line((d[p][0], sy * d[p][1]), (d[n][0], sy * d[n][1]), color))

    def draw_neurolucida():
        """  draw all neurolucida contours as SVG  """

        ### referential
        scene.add(Circle((0, 0), 5, (0, 0, 0), 1))
        # x
        scene.add(Line((0, 0), (20, 0), (0, 0, 0)))
        scene.add(Text((30, 0), "x", 16))
        # y
        scene.add(Line((0, 0), (0, sy * 20), (0, 0, 0)))
        scene.add(Text((0, sy * 30), "y", 16))

        ### draw contours
        draw_svg(contours['gran'], scene, colors['gran'])
        draw_svg(contours['glomint'], scene, colors['glom'])
        draw_svg(contours['glomext'], scene, colors['glom'])
        draw_svg(contours['line'], scene, (255, 0, 0))

    def center_angle():
        """    coordinate of center point and angle of line   """

        ### search intersection points between gran and line contours
        inter_points = []
        for p in contours['gran']:
            n = p + 1
            if n == len(contours['gran']):
               n = 0

            inter, inter_point = intersection(contours['gran'][p], contours['gran'][n], contours['line'][0], contours['line'][1])
            if DEBUG:
                print( inter, inter_point)

            if inter:
                inter_points.append(inter_point)
                if svg:
                    scene.add(Circle((inter_point[0], sy * inter_point[1]), 5, (255, 0, 0), 1))

        if svg: scene.write_svg()

        if verbose:
            print("\nIntersection points between granules contour and line: ", inter_points)
        if len(inter_points) != 2:
            print('\nIntersection points between granules ("gran" contour) and line must be 2!\n')
            return None, None, 2

        ### upper point of intersection between granules and line
        if inter_points[0][1] > inter_points[1][1]:   ### greater y
            xa, ya = inter_points[0][0], inter_points[0][1]
        else:
            xa, ya = inter_points[1][0], inter_points[1][1]

        ### coordinates of center of tissue
        xc = (inter_points[0][0] + inter_points[1][0]) / 2
        yc = (inter_points[0][1] + inter_points[1][1]) / 2
        if verbose: 
            print('Center of tissue:', xc, yc)

        if svg: scene.add(Circle((xc, sy * yc), 5, (0, 0, 255), 1))

        if svg:
            scene.write_svg()

        ### start angle (angle between line and vertical)
        #start_angle=asin(sqrt((xa-xc)**2)/sqrt((ya-yc)**2+(xa-xc)**2))

        start_angle = coord2angle((xc, yc), (xa, ya))

        if DEBUG:
            print('rough start angle', degrees(start_angle))

        if start_angle > pi / 2:
            start_angle -= 2 * pi

        if verbose:
            print('Start angle:', degrees(start_angle))

        return xc, yc, start_angle


    def sectors_contours():
        ### intersections
        inter_points = {'gran':{}, 'glomint':{}, 'glomext':{}}

        for sector in range(1, 360 // sector_degrees + 1):

            if verbose:
                print( 'Sector', sector)
            b = degrees(start_angle) + (sector - 1) * sector_degrees
            ### radius
            xm, ym = xc + sector_radius * sin(radians(b)), yc + sector_radius * cos(radians(b))
            if svg: scene.add(Line((xc, sy * yc), (xm, sy * ym), (128, 128, 128)))

            if svg: scene.write_svg()

            ### check intersections with contours
            for contour in ['gran', 'glomint', 'glomext']:

                for p in contours[contour]:
                    n = p + 1
                    if n == len(contours[contour]):
                       n = 0

                    """
                    if contour=='glomext' and sector==8:
                        print 'contour:',contour
                        print 'p:',p
                        print 'n:',n
                        print contours[contour][p]
                        print contours[contour][n]
                        print 'center',(xc,yc),'border',(xm,ym)
                    """

                    inter, inter_point = intersection(contours[contour][p], contours[contour][n], (xc, yc), (xm, ym))

                    """  
                    if contour=='glomext' and sector==8:
                        print 'inter:',inter
                        print inter_point
                    """

                    if inter:
                        inter_points[contour][sector] = inter_point
                        if verbose:
                            print('Intersection points:', inter_point)
                        if DEBUG and svg: scene.add(Circle((inter_point[0], sy * inter_point[1]), 5, (255, 0, 255), 1))  ### inter radius / gran

                    """
                    if contour=='glomext' and sector==8:
                        raw_input('key')   
                    """

        if verbose:
            for contour in ['gran', 'glomint', 'glomext']:
                print('Intersections between radius and "%s":' % contour, inter_points[contour], '\n')



        ### sectors contours
        c = 1
        area = {'gran':{}, 'glomint':{}, 'glomext':{}}

        for a in range(0, 360, sector_degrees):
            for contour in ['gran', 'glomint', 'glomext']:

                area[contour][c] = {}
                dummy = {}
                for p in contours[contour]:
                    angle = coord2angle((xc, yc), contours[contour][p]) - start_angle
                    if angle < 0:
                        angle += 2 * pi

                    if angle > 2 * pi:
                        angle -= 2 * pi

                    angle_deg = degrees(angle)
                    if angle_deg >= a and angle_deg < a + sector_degrees:
                        dummy[int(angle_deg)] = contours[contour][p]

                l = list(dummy.keys())  ### dummy for inserting points in angle ASC
                l.sort()
                #print(contour,c,l)

                for k in l:
                     area[contour][c][len(area[contour][c])] = dummy[k]

                n = c + 1
                if n == (360 // sector_degrees) + 1:
                     n = 1

                area[contour][c][len(area[contour][c])] = inter_points[contour][n]
                area[contour][c][len(area[contour][c])] = (xc, yc)
                area[contour][c][len(area[contour][c])] = inter_points[contour][c]

                if DEBUG and svg: draw_svg(area[contour][c], scene, (0, 255, 0))

            c += 1
        return area


    def markers_position():
        tot_mkr = {'glom':0, 'gran':0, 'epl':0}
        results = {}
        for b in range(0, 360, sector_degrees):
            results[len(results) + 1] = {'gran':0, 'glom':0, 'epl':0}
        ### markers angles
        for p in markers:

            ### inside granules?
            position = ''
            if inside_polygon(markers[p], contours['gran']):
                position += 'gran'

            ### inside glomerules?
            if not inside_polygon(markers[p], contours['glomint']) and inside_polygon(markers[p], contours['glomext']):
                position += 'glom'

            ### inside EPL
            if not inside_polygon(markers[p], contours['gran']) and inside_polygon(markers[p], contours['glomint']):
                position += 'epl'


            a = coord2angle((xc, yc), markers[p]) - start_angle
            if a < 0:
               a += 2 * pi

            sector = ((degrees(a) // sector_degrees) % (360 // sector_degrees)) + 1

            # draw marker

            if svg:
                 col = (0, 255, 0)
                 if position:
                     col = colors[position]
                 scene.add(Circle((markers[p][0], sy * markers[p][1]), 1, col, 1))
                 scene.add(Text((markers[p][0], sy * markers[p][1]), "%d s:%d" % (p + 1, sector), 5))

            if position:
                tot_mkr[position] += 1

            if position in results[sector]:
                results[sector][position] += 1
            else:
                results[sector][position] = 1

        return results, tot_mkr

    def write_file(outfile):
        """write a tsv data file with .txt extension"""

        if outfile:
            f = open(outfile, 'w')
        else:
            f = open(os.path.splitext(infile)[0] + '.txt', 'w')
            outfile = os.path.splitext(infile)[0]


        print( 'file:', os.path.splitext(os.path.basename(outfile))[0])


        #out='Region\tAnimal\tSlide\tx\tnb of markers\tx\tx\tx\tx\tx\tx\tx\tx\tx\tTotal\tSurface\n'
        out = 'Region\tAnimal\tLame\tCoupe\tGL\tEPL\tGrL\tRMS-OB\tONL\tCat6\tCat7\tCat8\tCat9\tCat10\tTotal\tSurface\t\t\t\t\t\t\t\t>>>\t42\n'

        out += 'Reg%(region)s\t\t\t%(coupe)s\t%(gran)d\t0\t0\t0\t0\t0\t0\t0\t0\t0\t%(tot_gran)d\t%(surface)d\n' % {'region':'1',
                                                                                                                 'coupe': os.path.splitext(os.path.basename(outfile))[0],
                                                                                                                 'gran': tot_mkr['gran'],
                                                                                                                 'tot_gran': tot_mkr['gran'],
                                                                                                                 'surface': polygon_area(contours['gran']),
                                                                                                                 'epl': tot_mkr['epl'],
                                                                                                                 'epl_surface': polygon_area(contours['glomint']) - polygon_area(contours['gran'])
                                                                                                                 }

        out += 'Reg%(region)s\t\t\t%(coupe)s\t%(glom)d\t0\t0\t0\t0\t0\t0\t0\t0\t0\t%(tot_glom)d\t%(surface)d\n' % {'region':'2',
                                                                                                                 'coupe': os.path.splitext(os.path.basename(outfile))[0],
                                                                                                                 'glom': tot_mkr['glom'],
                                                                                                                 'tot_glom': tot_mkr['glom'],
                                                                                                                 'surface': polygon_area(contours['glomext']) - polygon_area(contours['glomint'])}

        ct = 2
        for s in results:
            ct += 1
            out += 'Reg%(region)d\t\t\t%(coupe)s\t%(gran)d\t0\t0\t0\t0\t0\t0\t0\t0\t0\t%(tot_gran)d\t%(surface)d\n' % {'region':ct,
                                                                                            'coupe':os.path.splitext(os.path.basename(outfile))[0],
                                                                                            'gran': results[s]['gran'],
                                                                                            'tot_gran': results[s]['gran'],
                                                                                            'surface': polygon_area(aree['gran'][s]),
                                                                                            'epl': results[s]['epl'],
                                                                                           'epl_surface': polygon_area(aree['glomint'][s]) - polygon_area(aree['gran'][s])

                                                                                            }

        for s in results:
            ct += 1
            out += 'Reg%(region)d\t\t\t%(coupe)s\t%(glom)d\t0\t0\t0\t0\t0\t0\t0\t0\t0\t%(tot_glom)d\t%(surface)d\n' % {'region': ct,
                                                                                     'coupe': os.path.splitext(os.path.basename(outfile))[0],
                                                                                     'glom': results[s]['glom'],
                                                                                     'tot_glom': results[s]['glom'],
                                                                                     'surface': polygon_area(aree['glomext'][s]) - polygon_area(aree['glomint'][s])}

        if verbose:
            print(out)

        f.write(out)
        f.close()

        # EPL layer
        fepl = open(os.path.splitext(infile)[0] + '_epl.txt', 'w')
        out_epl = 'Region\tAnimal\tLame\tCoupe\tGL\tEPL\tGrL\tRMS-OB\tONL\tCat6\tCat7\tCat8\tCat9\tCat10\tTotal\tSurface\t\t\t\t\t\t\t\t>>>\t42\n'

        out_epl += 'Reg%(region)s\t\t\t%(coupe)s\t%(epl)d\t0\t0\t0\t0\t0\t0\t0\t0\t0\t%(epl)d\t%(epl_surface)d\n' % {'region': '1',
                                                                                                                 'coupe': os.path.splitext(os.path.basename(outfile))[0],
                                                                                                                 'gran': tot_mkr['gran'],
                                                                                                                 'tot_gran': tot_mkr['gran'],
                                                                                                                 'surface': polygon_area(contours['gran']),
                                                                                                                 'epl': tot_mkr['epl'],
                                                                                                                 'epl_surface': polygon_area(contours['glomint']) - polygon_area(contours['gran'])
                                                                                                                 }

        out_epl += 'Reg%(region)s\t\t\t%(coupe)s\t%(tot_glom)d\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t%(surface)d\n' % {'region': '2',
                                                                                           'coupe': os.path.splitext(os.path.basename(outfile))[0],
                                                                                           'glom': tot_mkr['glom'],
                                                                                           'tot_glom': tot_mkr['glom'],
                                                                                           'surface': polygon_area(contours['glomext']) - polygon_area(contours['glomint'])}

        # EPL sector count
        ct = 2
        for s in results:
            ct += 1
            out_epl += 'Reg%(region)d\t\t\t%(coupe)s\t%(epl)d\t0\t0\t0\t0\t0\t0\t0\t0\t0\t%(epl)d\t%(epl_surface)d\n' % {'region':ct,
                                                                                            'coupe':os.path.splitext(os.path.basename(outfile))[0],
                                                                                            'gran': results[s]['gran'],
                                                                                            'tot_gran': results[s]['gran'],
                                                                                            'surface': polygon_area(aree['gran'][s]),
                                                                                            'epl': results[s]['epl'],
                                                                                           'epl_surface': polygon_area(aree['glomint'][s]) - polygon_area(aree['gran'][s])

                                                                                            }

        ### glom as dummy value after EPL count
        for s in results:
            ct += 1
            out_epl += 'Reg%(region)d\t\t\t%(coupe)s\t%(glom)d\t0\t0\t0\t0\t0\t0\t0\t0\t0\t%(glom)d\t%(surface)d\n' % {'region':ct,
                                                                                     'coupe': os.path.splitext(os.path.basename(outfile))[0],
                                                                                     'glom': results[s]['glom'],
                                                                                     'tot_glom': results[s]['glom'],
                                                                                     'surface': polygon_area(aree['glomext'][s]) - polygon_area(aree['glomint'][s])
                                                                                     }


        if verbose:
            print(out_epl)


        fepl.write(out_epl)
        fepl.close()






    # check sector angle
    if 360 / sector_degrees != 360 // sector_degrees:
        return (5, '\nSector degrees must be a divisor of 360 not %d\n' % sector_degrees)

    contours, markers = parse_neurolucida_xml(infile)
    if contours == None:
        return (4, '\nNo contour found\n')

    if svg:
        if outfile:
            scene = Scene(os.path.splitext(outfile)[0])
        else:
            scene = Scene(os.path.splitext(infile)[0])

        draw_neurolucida()

        if svg: scene.write_svg()

    xc, yc, start_angle = center_angle()


    if xc == None and yc == None:
        return (2, 'Intersections between line and granules are <> 2')

    aree = sectors_contours()

    ### DEBUG
    #scene.write_svg()
    #sys.exit()

    results, tot_mkr = markers_position()

    if svg: scene.write_svg()

    if verbose:
        print('\n', results, '\n', tot_mkr)

    write_file(outfile)
    return (0, '')

if __name__ == '__main__':

    sector_degrees, verbose, svg, infile, outfile, DEBUG = check_options()

    if type(infile) == type([]):  ### file list
        for f in infile:
           n_error, text_error = main(sector_degrees, verbose, svg, f, outfile)
           if n_error:
               break
    else:
        n_error, text_error = main(sector_degrees, verbose, svg, infile, outfile)

    if n_error:
        print(text_error)
