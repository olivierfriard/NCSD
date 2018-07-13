#!/usr/bin/env python3

"""
Coordinate geometry module (http://www.mathopenref.com)

Copyright Olivier Friard 2010-2011
Universita' di Torino

This file is part of "Neurolucida Companion for Spatial Distribution Analysis".

"Neurolucida Companion for Spatial Distribution Analysis" is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or (at your option)
any later version.

"Neurolucida Companion for Spatial Distribution Analysis" is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with "Neurolucida Companion for Spatial Distribution Analysis"; see the file COPYING.TXT.  If not, write to
the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
"""


from __future__ import division

__version__ = '0.003'

DEBUG = False

from math import *


def DistancePointLine(P, M1, M2):
    """
    return:
    * distance from point P to segment M1-M2
    * coordinates of nearest point on segment
    """

    EPS = 0
    EPSEPS = EPS * EPS

    def SqLineMagnitude(x1, y1, x2, y2):
        return (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)

    px, py = P
    x1, y1 = M1
    x2, y2 = M2

    result = 0

    SqLineMag = SqLineMagnitude(x1, y1, x2, y2)
    if SqLineMag < EPSEPS:
        return - 1.0

    u = ((px - x1) * (x2 - x1) + (py - y1) * (y2 - y1)) / SqLineMag

    if (u < EPS) or (u > 1):
        ###  Closest point does not fall within the line segment,
        ###  take the shorter distance to an endpoint
        d1 = SqLineMagnitude(px, py, x1, y1)
        d2 = SqLineMagnitude(px, py, x2, y2)
        if d1 <= d2:
            result = d1
            ix, iy = x1, y1
        else:
            result = d2
            ix, iy = x2, y2

    else:

        ###  Intersecting point is on the line, use the formula
        ix = x1 + u * (x2 - x1)
        iy = y1 + u * (y2 - y1)
        result = SqLineMagnitude(px, py, ix, iy)


    return sqrt(result), (ix, iy)



def intersection(A, B, C, D):
    """
    line segments intersection with decimal precision
    return True and coordinates of intersection point otherwise
    return False,None
    """
    import decimal
    from decimal import getcontext
    getcontext().prec = 28

    Dec = decimal.Decimal
    xa, ya = Dec(str(A[0])), Dec(str(A[1]))
    xb, yb = Dec(str(B[0])), Dec(str(B[1]))
    xc, yc = Dec(str(C[0])), Dec(str(C[1]))
    xd, yd = Dec(str(D[0])), Dec(str(D[1]))

    """
    if round(yb)==-338 or round(yc)==-338:
        print 'segmenti',A,B,C,D
        print 'xa,ya,xb,yb',xa,ya,xb,yb
        print 'xc yc xd yd',xc,yc,xd,yd
    """

    ### check if first segment is vertical
    if xa == xb:
        slope = (yc - yd) / (xc - xd)
        intersept = yc - slope * xc
        #print 'slope, intersept',slope,intersept
        xm = xa
        ym = slope * xm + intersept

    ### check if second segment is vertical
    elif xc == xd:
        slope = (ya - yb) / (xa - xb)
        intersept = ya - slope * xa
        #print 'slope, intersept',slope,intersept
        xm = xc
        ym = slope * xm + intersept

    else:

        ### round Decimal result: .quantize(Dec('.001'), rounding=decimal.ROUND_DOWN) 
        xm = ((xd * xa * yc - xd * xb * yc - xd * xa * yb - xc * xa * yd + xc * xa * yb + xd * ya * xb + xc * xb * yd - xc * ya * xb) / (-yb * xd + yb * xc + ya * xd - ya * xc + xb * yd - xb * yc - xa * yd + xa * yc)).quantize(Dec('.001'), rounding=decimal.ROUND_DOWN)

        ym = ((yb * xc * yd - yb * yc * xd - ya * xc * yd + ya * yc * xd - xa * yb * yd + xa * yb * yc + ya * xb * yd - ya * xb * yc) / (-yb * xd + yb * xc + ya * xd - ya * xc + xb * yd - xb * yc - xa * yd + xa * yc)).quantize(Dec('.001'), rounding=decimal.ROUND_DOWN)


    xmin1, xmax1 = min(xa, xb), max(xa, xb)
    xmin2, xmax2 = min(xc, xd), max(xc, xd)

    ymin1, ymax1 = min(ya, yb), max(ya, yb)
    ymin2, ymax2 = min(yc, yd), max(yc, yd)


    if xm >= xmin1 and xm <= xmax1 and xm >= xmin2 and xm <= xmax2 and ym >= ymin1 and ym <= ymax1 and ym >= ymin2 and ym <= ymax2:
        return True, (float(xm), float(ym))
    else:
        return False, None


def polygon_area(poly):
    """
    area of polygon
    from http://www.mathopenref.com/coordpolygonarea.html
    """
    tot = 0
    for p in poly:
        x1, y1 = poly[p]
        n = (p + 1) % len(poly)
        x2, y2 = poly[n]
        tot += x1 * y2 - x2 * y1

    return abs(tot / 2)


def polygon_centroid(poly):
    """
    centroid of a polygon
    from http://local.wasp.uwa.edu.au/~pbourke/geometry/polyarea/
    """

    def centroid_polygon_area(poly):
        """
        area of polygon for centroid function
        from http://local.wasp.uwa.edu.au/~pbourke/geometry/polyarea
        
        return signed area!
        """
        tot = 0
        for p in poly:
            x1, y1 = poly[p]
            n = (p + 1) % len(poly)
            x2, y2 = poly[n]
            tot += x1 * y2 - x2 * y1

        return tot / 2

    cx = 0
    cy = 0
    for i in poly:

        j = (i + 1) % len(poly)
        p = poly[i][0] * poly[j][1] - poly[j][0] * poly[i][1]


        cx += (poly[i][0] + poly[j][0]) * p
        cy += (poly[i][1] + poly[j][1]) * p

    area = centroid_polygon_area(poly)

    cx = cx / (6 * area)
    cy = cy / (6 * area)

    return (cx, cy)


def inside_polygon(m, poly_d):
    """
    Determine if a point is inside a given polygon or not
    Polygon is a list of (x,y) pairs. This fuction
    returns True or False.  The algorithm is called
    "Ray Casting Method".
    """

    x, y = m

    poly = []
    for p in poly_d:
        poly.append(poly_d[p])

    n = len(poly)
    inside = False

    p1x, p1y = poly[0]
    for i in range(n + 1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside


def coord2angle2(m):
    a = atan2(m[1], m[0])
    return a


def coord2angle(ref, m):

    xc, yc = ref
    xk, yk = m
    sin_a = (xk - xc) / sqrt((xk - xc) ** 2 + (yk - yc) ** 2)
    cos_a = (yk - yc) / sqrt((xk - xc) ** 2 + (yk - yc) ** 2)

    if sin_a >= 0 and cos_a >= 0:
        a = asin(sin_a)
    if sin_a < 0 and cos_a >= 0:
        a = asin(sin_a) + 2 * pi
    if sin_a < 0 and cos_a < 0:
        a = -(asin(sin_a) - pi)
    if sin_a > 0 and cos_a < 0:
        a = acos(cos_a)

    return a

