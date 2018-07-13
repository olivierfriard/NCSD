#!/usr/bin/env python

"""
NCSD distance

a module for Neurolucida Companion for Spatial Distribution Analysis

Analyze markers in subcontours from arbitrary contour of a neurolucida XML file 

Copyright Olivier Friard 2010-2011
Universita' di Torino

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


version = '3'
"""
Module history

- 2011-05-12: fixed bug in draw_neurolucida function: added color control
- 2011-05-10: fixed bug in SVG drawing function   
"""


import sys
import os
from coordinate_geometry import *
import xml.etree.ElementTree as ET
from math import *
from svg import *
import glob
import glob

default_sector_degrees = 20
default_n_subcontours = 0

#contours_list = ['glomext', 'mitral', 'SVZ', 'line']

sy = -1
sector_radius = 2500
colors = {'line': (255, 0, 0), \
          'gran': (255, 228, 0), \
          'glom': (0, 255, 255), \
          'epl': (255, 0, 255), \
          'glomint': (0, 255, 255), \
          'glomext': (0, 255, 255), \
          'mitral': (0, 255, 0), \
          'SVZ': (0, 0, 255), \
          'subcontour': (200, 200, 200)          }


def check_options():
    """check command line options"""

    from optparse import OptionParser
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-a", "--angle", dest="angle", help="angle of sectors in degree (must be divisor of 360).\nThe default value is 20 degrees")
    parser.add_option("-c", "--contour", dest="ref_contour", help="reference contour")
    parser.add_option("-n", "--ncontours", dest="n_subcontours", help="number of subcontours")
    parser.add_option("-e", "--center_contour", dest="center_contour", help="cEnter contour for slide center determination")
    parser.add_option("-i", "--input", dest="input", help="Neurolucida XML file(s) (use \"*.xml\" to select more files)")
    parser.add_option("-o", "--output", dest="output", help="write results to file")
    parser.add_option("-v", "--verbose", action="store_true", default=False, dest="verbose", help="Print program messages")
    parser.add_option("-s", "--svg", action="store_true", default=False, dest="svg", help="Write SVG file")

    parser.add_option("-d", "--debug", action="store_true", default=False, dest="debug", help="Add information useful for debugging in output")

    (options, args) = parser.parse_args()

    if options.angle:
        sector_degrees = int(options.angle)
        if sector_degrees <= 0 or sector_degrees > 360:
            print('\nAngle must be >0 and <= 360 and a divisor of 360 and \n')
            sys.exit(19)

    else:
        sector_degrees = int(default_sector_degrees)

    if options.ref_contour:
        ref_contour = str(options.ref_contour)
    else:
        print( '\nNo contour selected!\nUse -h option for usage details\n')
        sys.exit(11)


    if options.n_subcontours:
        n_subcontours = int(options.n_subcontours)
    else:
        n_subcontours = int(default_n_subcontours)

    if options.center_contour:
        center_contour = str(options.center_contour)
    else:
        print( '\nNo contour selected for determination of slide center!\nUse -h option for usage details\n')
        sys.exit(15)


    verbose = options.verbose
    svg = options.svg
    DEBUG = options.debug

    if options.input:
        infile = options.input
        if '*' in infile:
            infile = glob.glob(infile)
    else:
        infile = input('Input XML file: ')
        if '*' in infile:
            infile = glob.glob(infile)

    if type(infile) == type(''):
        if not os.path.exists(infile):
            print( '\n"%s" not found!\nUse -h option for usage details\n' % infile)
            sys.exit(1)

    if options.output:
        outfile = options.output
    else:
        outfile = ''
    return (sector_degrees, ref_contour, center_contour, n_subcontours, verbose, svg, infile, outfile, DEBUG)



def main(sector_degrees, ref_contour, center_contour, n_subcontours, verbose, svg, infile, outfile):

    def parse_neurolucida_xml(infile):
        """
        parse neurolucida XML file and returns list of contours and list of markers
        """
        ### create dictionary
        contours = {}
        markers = {}

        tree = ET.parse(infile)
        root = tree.getroot()

        for contour in root.getiterator('contour'):

            name = contour.findtext('name')

            if verbose:
                print('found "%s" contour' % name)
            if not name in contours:
                contours[name] = {}

            for point in contour.getiterator('point'):
                if name in contours:
                    contours[name][len(contours[name])] = (float(point.get('x')), float(point.get('y')))

        for c in contours:
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

    def draw_svg(d, scene, color, draw_number=True , draw_dot=True):

        ### mark 1st point of contour
        if DEBUG:
            scene.add(Circle((d[0][0], sy * d[0][1]), 5, color, 1))
            ### point number
            #scene.add(Text((d[0][0]+4,sy*d[0][1]),"%d" % (1),10))

        for p in d:
            n = p + 1
            if n == len(d):
               n = 0

            ### draw point of contour
            if DEBUG:
               if draw_dot:
                   scene.add(Circle((d[p][0], sy * d[p][1]), 2, (0, 0, 0), 1))

               ### point number
               if draw_number:
                   scene.add(Text((d[p][0] + 4, sy * d[p][1]), "%d" % (p + 1), 10))

            scene.add(Line((d[p][0], sy * d[p][1]), (d[n][0], sy * d[n][1]), color))


    def draw_neurolucida(contours_list):
        """
        draw all neurolucida contours as SVG
        """

        ### referential
        scene.add(Circle((0, 0), 5, (0, 0, 0), 1))
        # x
        scene.add(Line((0, 0), (20, 0), (0, 0, 0)))
        scene.add(Text((30, 0), "x", 16))
        # y
        scene.add(Line((0, 0), (0, sy * 20), (0, 0, 0)))
        scene.add(Text((0, sy * 30), "y", 16))

        ### draw contours

        for c in contours_list:

            ### check if contour has defined color
            if c in colors:
                col = colors[c]
            else:
                col = (0, 255, 255)

            draw_svg(contours[c], scene, col)

    def center_angle(center_contour):
        """
        coordinate of center point and angle of line
        """

        ### search intersection points between gran and line contours
        #intersection_contour = 'SVZ'
        #intersection_contour = ref_contour


        inter_points = []
        for p in contours[center_contour]:
            n = p + 1
            if n == len(contours[center_contour]):
               n = 0

            if DEBUG:
                print( 'p n:', p, n)
                print( contours[center_contour][p], contours[center_contour][n])
                print( 'line', contours['line'][0], contours['line'][1])

            inter, inter_point = intersection(contours[center_contour][p], contours[center_contour][n], contours['line'][0], contours['line'][1])
            if DEBUG:
                print(inter, inter_point)

            if inter:
                inter_points.append(inter_point)
                if svg:
                    scene.add(Circle((inter_point[0], sy * inter_point[1]), 5, (255, 0, 0), 1))

        if svg: scene.write_svg()

        if verbose:
            print( "\nIntersection points between granules contour and line: ", inter_points)
        if len(inter_points) != 2:
            print('\nIntersection points between "' + center_contour + '" contour and line must be 2!\n')
            return None, None, 2

        ### upper point of intersection between contour and line
        if inter_points[0][1] > inter_points[1][1]:   ### greater y
            xa, ya = inter_points[0][0], inter_points[0][1]
        else:
            xa, ya = inter_points[1][0], inter_points[1][1]

        ### coordinates of center of tissue
        xc = (inter_points[0][0] + inter_points[1][0]) / 2
        yc = (inter_points[0][1] + inter_points[1][1]) / 2
        if verbose:
            print('Center of tissue:', xc, yc)

        if svg: scene.add(Circle((xc , sy * yc), 5, (0, 0, 0), 1))

        if svg: scene.write_svg()

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


    def sectors_contours(ref_contour, contours):

        ### intersections
        inter_points = {}
        #inter_points [ref_contour] = {}

        for contour in contours:
            if ref_contour in contour:
                inter_points [contour] = {}

        for sector in range(1, 360 // sector_degrees + 1):

            if verbose:
                print( 'Sector', sector)
            b = degrees(start_angle) + (sector - 1) * sector_degrees
            ### radius
            xm, ym = xc + sector_radius * sin(radians(b)), yc + sector_radius * cos(radians(b))
            if svg: scene.add(Line((xc, sy * yc), (xm, sy * ym), (128, 128, 128)))

            if svg: scene.write_svg()

            ### check intersections with contours
            for contour in contours:
                if not ref_contour in contour:
                    continue

                for p in contours[contour]:
                    n = p + 1
                    if n == len(contours[contour]):
                       n = 0

                    inter, inter_point = intersection(contours[contour][p], contours[contour][n], (xc, yc), (xm, ym))

                    if inter:
                        inter_points[contour][sector] = inter_point
                        if verbose:
                            print('Intersection points:', inter_point)
                        if DEBUG and svg: scene.add(Circle((inter_point[0], sy * inter_point[1]), 5, (255, 0, 255), 1))  ### inter radius / gran

        if verbose:
            for contour in contours:
                if ref_contour in contour :
                    print('Intersections between radius and "%s":' % contour, inter_points[contour], '\n')


        ### sectors contours
        c = 1

        ### init dictionary for area of contour in each sector
        area = {}
        for contour in contours:
            if ref_contour in contour:
                area [contour] = {}

        for a in range(0, 360, sector_degrees):

            for contour in contours:

                if not ref_contour in contour:
                    continue

                if DEBUG:
                    print('contour:', contour)

                area[contour][c] = {}
                dummy = {}
                for p in contours[contour]:
                    if DEBUG:
                        print( 'p:', p + 1)
                    angle = coord2angle((xc, yc), contours[contour][p]) - start_angle
                    if DEBUG: 
                        print('angle:', degrees(angle))
                    if angle < 0:
                        angle += 2 * pi

                    if angle > 2 * pi:
                        angle -= 2 * pi

                    if DEBUG:
                        print('angle modif:', degrees(angle))
                        print( 'a a+sector:', a, a + sector_degrees)

                    angle_deg = degrees(angle)
                    if angle_deg >= a and angle_deg < a + sector_degrees:
                        dummy[int(angle_deg)] = contours[contour][p]

                    #raw_input('key')    

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


    def markers_position(ref_contour, n_subcontours, contours):
        ''' count markers in ref_contour and subcontour of ref_contour'''

        def get_sector(marker):
            ''' determine sector of marker'''
            a = coord2angle((xc, yc), marker) - start_angle
            if a < 0:
               a += 2 * pi

            if DEBUG:
                print('a', degrees(a))

            return int(((degrees(a) // sector_degrees) % (360 // sector_degrees)) + 1)


        # init dictionary for all results
        results = {}


        for sector in range(1, 360 // sector_degrees + 1):

            results[sector] = {}

            for contour in contours:
                if ref_contour in contour:
                    results[sector] [contour] = 0
            """
            for i in range(1, n_subcontours):
                results[sector] [ref_contour + '|' + str(i)] = 0
            """


        # init dictionary for number of markers
        tot_mkr = {}
        for contour in contours:
            if ref_contour in contour:
                tot_mkr[contour] = 0

        """
        for i in range(1, n_subcontours + 1):
            tot_mkr[ref_contour + '|' + str(i)] = 0
        """

        # markers angles
        for p in markers:

            position = ''
            # check if marker in ref_contour

            if inside_polygon(markers[p], contours[ref_contour]):
                   tot_mkr[ref_contour] += 1
                   results[get_sector(markers[p])][ref_contour] += 1

            # check if marker in subcontour
            for i in range(1, n_subcontours + 1):
                if i > 1:
                    if inside_polygon(markers[p], contours[ref_contour + '|' + str(i)]) and not inside_polygon(markers[p], contours[ref_contour + '|' + str(i - 1)]):
                        #position += ref_contour + '|' + str(i)
                        tot_mkr[ref_contour + '|' + str(i)] += 1
                        results[get_sector(markers[p])][ref_contour + '|' + str(i)] += 1

                else:
                    if inside_polygon(markers[p], contours[ref_contour + '|' + str(i)]):
                        #position += ref_contour + '|' + str(i)
                        tot_mkr[ref_contour + '|' + str(i)] += 1
                        results[get_sector(markers[p])][ref_contour + '|' + str(i)] += 1

            ### sector?

            if DEBUG:
                print('coord angle', degrees(coord2angle((xc, yc), markers[p])))
                print('start_angle', degrees(start_angle))
            sector = get_sector(markers[p])

            if verbose:
                print('Marker: %d   sector: %d  in %s' % (p + 1, get_sector(markers[p]), position))    #degrees(acos(cos_a)) )

            ### draw marker

            if svg:
                 col = (0, 255, 0)
                 scene.add(Circle((markers[p][0], sy * markers[p][1]), 1, col, 1))
                 scene.add(Text((markers[p][0], sy * markers[p][1]), "%d s:%d" % (p + 1, get_sector(markers[p])), 5))

        return results, tot_mkr

    def write_file(ref_contour, n_subcontours, outfile):
        """write a tsv data file with .txt extension"""

        if os.path.isdir(outfile):
            print( outfile + ' is a directory')
            return False
        if outfile:
            f = open(outfile, 'w')
        else:
            f = open(os.path.splitext(infile)[0] + '.txt', 'w')
            outfile = os.path.splitext(infile)[0]

        if DEBUG:
            print('file:', os.path.splitext(os.path.basename(outfile))[0])

        out = 'Region\tTotal markers\t'
        for s in results:
            out += 'Markers in sector ' + str(s) + '\t'
        for s in results:
            out += 'Area of sector ' + str(s) + '\t'

        out = out[0:-1] + '\n'


        out += ref_contour + '\t'   #+ '\t\t\t\t\t\t'



        # number of markers for ref contour
        #total
        out += str(tot_mkr[ref_contour]) + '\t'

        ### by sector
        for s in results:
            out += str(results[s][ref_contour]) + '\t'

        ### sector area for reference contour
        for s in results:
            out += str(int(polygon_area(aree[ref_contour][s]))) + '\t'

        out = out[0:-1] + '\n'

        for i in range(1, n_subcontours + 1):
             out += '%(ref_contour)s %(index)d\t%(tot_mkr)d\t' % {'ref_contour':ref_contour, 'index':i, 'tot_mkr':tot_mkr[ref_contour + '|' + str(i)]}

             ### markers number by sector
             for s in results:
                 out += str(results[s][ref_contour + '|' + str(i)]) + '\t'

             ### area subcontour by sector
             for s in results:
                 if i > 1:
                     out += str(int(polygon_area(aree[ref_contour + '|' + str(i)][s]) - polygon_area(aree[ref_contour + '|' + str(i - 1)][s]))) + '\t'
                 else:
                     out += str(int(polygon_area(aree[ref_contour + '|' + str(i)][s]))) + '\t'


             out = out[0:-1] + '\n'

        print(out)

        f.write(out)
        f.close()


        return True


    """main body"""
    ### check sector angle
    if 360 / sector_degrees != 360 // sector_degrees:
        return (5, '\nSector degrees must be a divisor of 360 not %d\n' % sector_degrees)

    contours, markers = parse_neurolucida_xml(infile)

    if contours == None:
         return (4, '\nNo contour found\n')

    if not ref_contour in contours:

        return (12, '\n%s contour not found in Neurolucida XML file!\n' % (ref_contour) + 'Existing contours are:\n' + '\n'.join(contours.keys()) + '\n')

    if not center_contour in contours:

        return (15, '\n%s contour not found in Neurolucida XML file!\n' % (center_contour) + 'Existing contours are:\n' + '\n'.join(contours.keys()) + '\n')


    if len(contours[ref_contour]) < 3:
        return (13, '\nThe reference contour has not enougth points\n')

    if svg:
        if outfile:
            scene = Scene(os.path.splitext(outfile)[0])
        else:
            scene = Scene(os.path.splitext(infile)[0])

        draw_neurolucida(contours)

        scene.write_svg()

    xc, yc, start_angle = center_angle(center_contour)

    if xc == None and yc == None:
        return (2, 'Intersections between line and center contour are != 2')

    ### contour design

    ### init sub contours

    for i in range(1, n_subcontours + 1):
        contours[ref_contour + '|' + str(i)] = {}

    for p in contours[ref_contour]:
            n = p + 1
            if n == len(contours[ref_contour]):
               n = 0

            for i in range(1, n_subcontours + 1):
                contours[ref_contour + '|' + str(i)][p] = (xc + (contours[ref_contour][p][0] - xc) * (i / n_subcontours), \
                                                           yc + (contours[ref_contour][p][1] - yc) * (i / n_subcontours))

    if svg:
        for i in range(1, n_subcontours + 1):
            draw_svg(contours[ref_contour + '|' + str(i)], scene, colors['subcontour'], False, False) ### no number, no dot

    aree = sectors_contours(ref_contour, contours)

    if DEBUG:
        print( 'Aree:\n', aree)


    results, tot_mkr = markers_position(ref_contour, n_subcontours, contours)

    if DEBUG:
        print( 'results', results)
        print( 'tot_mkr', tot_mkr)

    if svg: scene.write_svg()

    if write_file(ref_contour, n_subcontours, outfile):
        return (0, '')
    else:
        return (20, 'Error writing result file')


### MAIN
if __name__ == '__main__':

    sector_degrees, ref_contour, center_contour, n_subcontours, verbose, svg, infile, outfile, DEBUG = check_options()

    if type(infile) == type([]):
        for f in infile:
           n_error, text_error = main(sector_degrees, ref_contour, center_contour, n_subcontours, verbose, svg, f, outfile)
           if n_error:
               break
    else:
        n_error, text_error = main(sector_degrees, ref_contour, center_contour, n_subcontours, verbose, svg, infile, outfile)

    if n_error:
        print(text_error)
