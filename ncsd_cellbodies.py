#!/usr/bin/env python

"""
NCSD cell bodies

a module for Neurolucida Companion for Spatial Distribution Analysis

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

from __future__ import division

version = '4'

import sys, os
from coordinate_geometry import *
import xml.etree.ElementTree as ET
from math import *
from svg import *
import os
import glob

default_sector_degrees = 20

sy = -1
sector_radius = 2500
colors = {'gran':(255, 228, 0), 'glom':(0, 255, 255), 'cell':(100, 100, 100)}

def check_options():
    """check command line options"""
    from optparse import OptionParser
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-a", "--angle", dest="angle", help="angle of sectors in degree (must be divisor of 360).\nThe default value is 20 degrees")
    parser.add_option("-c", "--contour", dest="ref_contour", help="reference contour")
    parser.add_option("-e", "--center_contour", dest="center_contour", help="cEnter contour for slide center determination")
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

    if options.ref_contour:
        ref_contour = str(options.ref_contour)
    else:
        print('\nNo contour selected!\nUse -h option for usage details\n')
        sys.exit(11)

    if options.center_contour:
        center_contour = str(options.center_contour)
    else:
        print('\nNo contour selected for determination of slide center!\nUse -h option for usage details\n')
        sys.exit(15)


    verbose = options.verbose
    svg = options.svg
    DEBUG = options.debug

    if options.input:
        infile = options.input
        if '*' in infile:
            infile = glob.glob(infile)
    else:
        infile = raw_input('Input XML file: ')
        if '*' in infile:
            infile = glob.glob(infile)

    if type(infile) == type(''):
        if not os.path.exists(infile):
            print('\n"%s" not found!\n' % infile)
            sys.exit(1)

    if options.output:
        outfile = options.output
    else:
        outfile = ''
    return (sector_degrees, ref_contour, center_contour, verbose, svg, infile, outfile, DEBUG)
#return (sector_degrees, ref_contour, n_subcontours, verbose, svg, infile, outfile, DEBUG)



def main(sector_degrees, ref_contour, center_contour, verbose, svg, infile, outfile):

    def parse_neurolucida_xml(infile):
        """
        parse neurolucida XML file and returns list of contours and list of markers
        """

        #contours = {'glomint':{}, 'glomext':{}, 'gran':{}, 'line':{}}
        ### create dictionary
        contours = {}
        #markers = {}


        tree = ET.parse(infile)
        root = tree.getroot()

        for contour in root.getiterator('contour'):

            name = contour.findtext('name')

            if verbose:
                print('found "%s" contour' % name)

            if not name in contours.keys():
                contours[name] = {}

            for point in contour.getiterator('point'):
                contours[name][len(contours[name])] = (float(point.get('x')), float(point.get('y')))

        for c in contours:
            if verbose:
                print('"%s" contour vertices number: %d' % (c, len(contours[c])))
            if not len(contours[c]):
                print('\n"%s" contour is invalid (%d vertice(s))\n' % (c, len(contours[c])))
                return None, None

        return contours


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
        for c in contours_list:

            ### check if contour has defined color
            if c in colors:
                col = colors[c]
            else:
                col = (0, 255, 255)

            draw_svg(contours[c], scene, col)

    def center_angle(center_contour):
        """    coordinate of center point and angle of line   """

        ### search intersection points between gran and line contours
        inter_points = []
        for p in contours[center_contour]:
            n = p + 1
            if n == len(contours[center_contour]):
               n = 0

            if DEBUG:
                print('p n:', p, n)
                print(contours[center_contour][p], contours[center_contour][n])
                print('line', contours['line'][0], contours['line'][1])

            inter, inter_point = intersection(contours[center_contour][p], contours[center_contour][n], contours['line'][0], contours['line'][1])
            if DEBUG:
                print(inter, inter_point)

            if inter:
                inter_points.append(inter_point)
                if svg: scene.add(Circle((inter_point[0], sy * inter_point[1]), 5, (255, 0, 0), 1))

        if svg: scene.write_svg()

        if verbose:
            print("\nIntersection points between " + center_contour + " contour and line: ", inter_points)
        if len(inter_points) != 2:
            print('\nIntersection points between "' + center_contour + '" contour and line must be 2!\n')
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

    def draw_sectors():
        ### intersections

        for sector in range(1, 360 // sector_degrees + 1):

            if DEBUG:
                 print('Sector', sector)
            b = degrees(start_angle) + (sector - 1) * sector_degrees
            ### radius
            xm, ym = xc + sector_radius * sin(radians(b)), yc + sector_radius * cos(radians(b))
            scene.add(Line((xc, sy * yc), (xm, sy * ym), (128, 128, 128)))

        return



    def write_file(outfile):
        """write a tsv data file"""

        out = 'Cell ID\tsector\tArea\tDistance from ' + ref_contour + '\n'

        centroids_keys = centroids.keys()
        centroids_keys.sort()
        for centroid in centroids_keys:
            out += '%(centroid)s\t%(sect)d\t%(area).2f\t%(dist).2f\n' % {'centroid': centroid, 'sect': centroids_sector[centroid], 'area': centroids_area[centroid], 'dist': distances[centroid][0]}
            #centroid + '\t' + str(int(centroids_sector[centroid])) + '\t' + str(centroids_area[centroid]) + '\t' + str(distances[centroid][0]) + '\n'
            ### %.1f' % round(n, 1) gives you '5.6'
        if verbose:
            print(out)

        if outfile:
            f = open(outfile, 'w')
        else:
            f = open(os.path.splitext(infile)[0] + '.tsv', 'w')
        f.write(out)
        f.close()
        return True


    """main function"""
    ### check sector angle
    if 360 / sector_degrees != 360 // sector_degrees:
        print('Sector degrees must be a divisor of 360 not %d' % sector_degrees)
        return 5

    ### parse xml file for contours and markers

    contours = parse_neurolucida_xml(infile)

    if contours == None:
        return (4, '\nNo contour found\n')

    if not ref_contour in contours:

        return (12, '\n%s contour not found in Neurolucida XML file!\n' % (ref_contour) + 'Existing contours are:\n' + '\n'.join(contours.keys()) + '\n')

    if not center_contour in contours:

        return (15, '\n%s contour not found in Neurolucida XML file!\n' % (center_contour) + 'Existing contours are:\n' + '\n'.join(contours.keys()) + '\n')


    if svg:
        if outfile:
            scene = Scene(os.path.splitext(outfile)[0])
        else:
            scene = Scene(os.path.splitext(infile)[0])

        draw_neurolucida(contours)

        scene.write_svg()

    # determine centroids of cells contours
    centroids = {}
    for c in contours.keys():
        if 'cell' in c:
            centroids[c] = polygon_centroid(contours[c])

    # aree of cells
    centroids_area = {}
    for c in contours.keys():
        if 'cell' in c:
            centroids_area[c] = polygon_area(contours[c])
    if verbose:
        print('\nCells area', centroids_area)

    # determine distance between glomint contour and centroid's cells
    distances = {}
    for centroid in centroids:
        mind = 1e9
        minip = (1e9, 1e9)
        if DEBUG:
            print('centroid', centroid)
        for p in contours[ref_contour]:

            n = (p + 1) % len(contours[ref_contour])

            d, ip = DistancePointLine(centroids[centroid], contours[ref_contour][p] , contours[ref_contour][n])

            if d < mind:
                mind = d
                minip = ip

        distances[centroid] = (mind, minip)

    if verbose:
        print('\nDistances from ' + ref_contour + ' contour', distances)

    # draw distances between centroid and ref_contour 
    if svg:
        for centroid in centroids:
            xc, yc = centroids[centroid]
            xd, yd = distances[centroid][1]
            scene.add(Line((xc, sy * yc), (xd, sy * yd), (154, 154, 154)))


    # determine slide center and angle 
    xc, yc, start_angle = center_angle(center_contour)

    if xc == None and yc == None:
        return 2 # error 2: intersections between line and granules are <> 2

    # draw sectors
    if svg:
        draw_sectors()


    # centroids sector
    centroids_sector = {}
    for centroid in centroids:
        if DEBUG:
            print('centroid', centroid)
        a = coord2angle((xc, yc), centroids[centroid]) - start_angle
        if a < 0:
            a += 2 * pi

        if DEBUG:
            print('centroid angle', centroid, degrees(a))
        centroids_sector[centroid] = ((degrees(a) // sector_degrees) % (360 // sector_degrees)) + 1
        if DEBUG:
            print('centroids sector', centroid, centroids_sector[centroid])
        if svg:
            scene.add(Circle((centroids[centroid][0], sy * centroids[centroid][1]), 2, (255, 0, 0), 1))
            scene.add(Text((centroids[centroid][0], sy * centroids[centroid][1]), centroid + ' ' + str(centroids_sector[centroid]), 12))

    if verbose:
        print('\nCentroids_sector', centroids_sector)
    #raw_input('key')    

    if svg: scene.write_svg()

    write_file(outfile)

    return (0, '') # no error

if __name__ == '__main__':

    sector_degrees, ref_contour, center_contour, verbose, svg, infile, outfile, DEBUG = check_options()

    if type(infile) == type([]):
        for f in infile:
           n_error, text_error = main(sector_degrees, ref_contour, center_contour, verbose, svg, f, outfile)
           if n_error:
               break
    else:
        n_error, text_error = main(sector_degrees, ref_contour, center_contour, verbose, svg, infile, outfile)

    if n_error:
        print(text_error)
