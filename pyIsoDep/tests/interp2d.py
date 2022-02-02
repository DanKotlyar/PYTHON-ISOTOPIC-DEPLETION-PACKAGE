# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 09:57:30 2022

@author: dkotlyar6
"""

import bisect
import numpy


def _interp2D(x, xvalues, y, yvalues, Z):

    if x < min(xvalues) or x > max(xvalues):
        raise ValueError(
            "x must be between {} and {}, not {}".format(
                min(xvalues), max(xvalues), x
            )
        )
    if y < min(yvalues) or y > max(yvalues):
        raise ValueError(
            "y must be between {} and {}, not {}".format(
                min(yvalues), max(yvalues), y
            )
        )

    # Find the extreme cases (P,T)min and (P,T)max
    idx00 = numpy.intersect1d(numpy.where(xvalues <= x), numpy.where(
            yvalues <= y), return_indices=False)[-1]
    idx11 = numpy.intersect1d(numpy.where(xvalues >= x), numpy.where(
            yvalues >= y), return_indices=False)[0]

    # (P,T) exist and there is no need to interpolate
    if idx00 == idx11:
        return Z[idx00]
    # same P[MPa], but different T[K]
    if xvalues[idx00] == xvalues[idx11]:
        ypts = yvalues[idx00], yvalues[idx11]
        zpts = Z[idx00], Z[idx11]
        return _local1DInterp(y, ypts, zpts)
    # same T[K], but different P[MPa]
    elif yvalues[idx00] == yvalues[idx11]:
        xpts = xvalues[idx00], xvalues[idx11]
        zpts = Z[idx00], Z[idx11]
        return _local1DInterp(x, xpts, zpts)

    zvalues = [
        [Z[idx00], Z[idx00+1]],
        [Z[idx11-1], Z[idx11]],
    ]

    xpts = xvalues[idx00], xvalues[idx11]
    ypts = yvalues[idx00], yvalues[idx11]

    return _bilinear2D(
        x,
        y,
        xpts,
        ypts,
        zvalues,
    )


def _bilinear2D(x, y, xv, yv, zm):
    denom = (xv[1] - xv[0]) * (yv[1] - yv[0])
    xlead = [xv[1] - x, x - xv[0]]
    ytail = [yv[1] - y, y - yv[0]]
    prod = numpy.matmul(zm, ytail)
    return numpy.matmul(xlead, prod) / denom


def _interp1D(x, xvalues, xdesc, yvalues):
    if x < min(xvalues) or x > max(xvalues):
        raise ValueError(
            "{} must be between {} and {}, not {}".format(
                xdesc, min(xvalues), max(xvalues), x
            )
        )
    # Find index that is closest to requested value
    index = bisect.bisect_left(xvalues, x)
    if xvalues[index] == x:
        return yvalues[index]
    return _local1DInterp(x, xvalues[index:index+2],
                          yvalues[index:index+2])


def _local1DInterp(c, x, y):
    assert len(x) == len(y)
    slope = (y[1] - y[0]) / (x[1] - x[0])
    return y[0] + slope * (c - x[0])
