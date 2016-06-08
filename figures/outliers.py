from __future__ import division

import sys
from numpy import abs, floor, array, std, mean, nan
from collections import defaultdict

from rpy2.robjects import DataFrame, FloatVector, StrVector, FactorVector, r
from rpy2.robjects.packages import importr
from handy_r import quantile
quantreg = importr('quantreg')

__author__ = 'Ben Kaehler'
__copyright__ = 'Copyright 2014, Ben Kaehler'
__credits__ = ['Ben Kaehler']
__license__ = 'GPL'
__maintainer__ = 'Ben Kaehler'
__email__ = 'benjamin.kaehler@anu.edu.au'
__status__ = 'Production'
__version__ = '0.0.4-dev'

def qlim(x, y):
    df = DataFrame({'x':FloatVector(x), 'y':FloatVector(y)})
    rq = quantreg.rq('y ~ x', df, tau=FloatVector((0.25, 0.5, 0.75)))
    print rq.rx2('coefficients')
    fv = array(rq.rx2('fitted.values'))
    return min(fv[:,1]) - 2*max(fv[:,1] - fv[:,0]), \
            2*max(fv[:,2] - fv[:,1]) + max(fv[:,1])

def qcrop(xlist, ylist, labels=None):
    if labels is None:
        labels = map(str, range(len(xlist)))
    x = []
    y = []
    xcrop = []
    ycrop = []
    facet = []
    for i, (onex, oney) in enumerate(zip(xlist, ylist)):
        ymin, ymax = qlim(onex, oney)
        cropx, cropy = zip(*[(nan, nan) if vy > ymax or vy < ymin else (vx, vy)
                for vx, vy in zip(onex, oney)])
        xcrop += cropx
        ycrop += cropy
        x += onex
        y += oney
        facet += [labels[i]]*len(onex)

    df = DataFrame({'x': FloatVector(x), 'y': FloatVector(y),
            'xcrop': FloatVector(xcrop) , 'ycrop': FloatVector(ycrop),
            'facet': FactorVector(StrVector(facet), levels=StrVector(labels))})
    return df

def qlim1(x, nq=4.):
    qs = array(r.quantile(FloatVector(x), probs=FloatVector((0.25,0.5,0.75))))
    return qs[1] - nq*(qs[1]-qs[0]), qs[1] + nq*(qs[2]-qs[1])

def qcrop2(xlist, ylist, labels=None, nq=4.):
    if labels is None:
        labels = map(str, range(len(xlist)))
    x = []
    y = []
    xcrop = []
    ycrop = []
    facet = []
    for i, (onex, oney) in enumerate(zip(xlist, ylist)):
        xmin, xmax = qlim1(onex, nq)
        ymin, ymax = qlim1(oney, nq)
        cropx, cropy = zip(*[(nan, nan) 
            if vy > ymax or vy < ymin or vx < xmin or vx > xmax else (vx, vy)
                for vx, vy in zip(onex, oney)])
        xcrop += cropx
        ycrop += cropy
        x += onex
        y += oney
        facet += [labels[i]]*len(onex)

    df = DataFrame({'x': FloatVector(x), 'y': FloatVector(y),
            'xcrop': FloatVector(xcrop) , 'ycrop': FloatVector(ycrop),
            'facet': FactorVector(StrVector(facet), levels=StrVector(labels))})
    return df

def combine_crops(crop1, crop2):
    if crop1 == (None, None):
        return crop2
    if crop2 == (None, None):
        return crop1
    return min(crop1[0],crop2[0]), max(crop1[1], crop2[1])

def qcrop_facet_grid(points, nq=4.):
    xlabels = list(set([row[2] for row in points]))
    ylabels = list(set([row[3] for row in points]))
    subsets = defaultdict(list)
    for row in points:
        subsets[row[2:]].append(row[:2])
    xcrops = defaultdict(lambda: (None, None))
    for xlabel in xlabels:
        for label, sample in subsets.items():
            if not xlabel in label:
                continue
            xsample = [p[0] for p in sample]
            xcrop = qlim1(xsample, nq)
            xcrops[xlabel] = combine_crops(xcrops[xlabel], xcrop)
    ycrops = defaultdict(lambda: (None, None))
    for ylabel in ylabels:
        for label, sample in subsets.items():
            if not ylabel in label:
                continue
            ysample = [p[1] for p in sample]
            ycrop = qlim1(ysample, nq)
            ycrops[ylabel] = combine_crops(ycrops[ylabel], ycrop)
    cropped_points = []
    for row in points:
        xcrop = xcrops[row[2]]
        ycrop = ycrops[row[3]]
        if row[0] > xcrop[0] and row[0] < xcrop[1] and \
                row[1] > ycrop[0] and row[1] < ycrop[1]:
            cropped_points.append(row + row[:2])
        else:
            cropped_points.append(row + (nan, nan))
    return cropped_points

def main():
    pass

if __name__ == '__main__':
    sys.exit(main())
