import math
import numpy
import json
import sys
from pprint import pprint
import pyproj
import itertools
from vincenty import vincenty
import operator


# Test data; these could be your inputs. This is actually scraped of a website
# and not generated. The website in question only supplies real positions for
# a subset of search results but we use those to validate accuracy.
real = [48.14696,11.73628]
coords = [11.580218296782743,48.13918807515703,11.477],[11.517138039698214,48.18768612414934,16.672],[11.5608050377351,48.153738136990604,12.889],[11.617346446831283,48.11674732163194,9.331],[11.699343887674134,48.049757497992914,10.99],[11.552912258276171,48.153474425856174,13.465],[11.53954333934099,48.07802969169499,16.294],[11.697149289797098,48.09768441519491,6.115],[11.552367644497043,48.090519951888496,14.844],[11.590195038225199,48.088195668989414,12.504],[11.508514050229243,48.1922851309257,17.418],[11.478152945885506,48.1040735748235,19.512],[11.702692394012457,48.20727395298781,7.053],[11.63383506324293,48.13446046433374,7.637],[11.518997884978273,48.22400790883851,18.022],[11.597246216512007,48.0816412935017,12.462],[11.663962187263063,48.07012805230099,9.952],[11.715425377826055,48.12508875839166,2.842],[11.694984339584895,48.05326530794745,10.705],[11.684659307026715,48.094214333519645,6.91],[11.606752133062525,48.166948860628295,9.746],[11.703320847275219,48.19492282616845,5.784],[11.655273933378785,48.12798187910728,6.295],[11.568127687109664,48.07484573419233,14.653],[11.517470445042838,48.156072308432634,16.076],[11.548151207567555,48.0824097832868,15.511],[11.519083704698566,48.158870724575216,15.979],[11.60935979193895,48.22311257449991,12.495],[11.528556605113065,48.07641081289575,17.092],[11.478680391889089,48.12322489032875,19.073],[11.604630589559145,48.10073011293635,10.906]

# This is used to correct distances provided by the scraped data and was
# determined experimentally. You might not need or want this.
MagicDistanceFactor = 1.0142

# Projections used in the code.
ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')

def trilaterate(points):
    LatA = points[0][1]
    LonA = points[0][0]
    DistA = points[0][2]*MagicDistanceFactor
    LatB = points[1][1]
    LonB = points[1][0]
    DistB = points[1][2]*MagicDistanceFactor
    LatC = points[2][1]
    LonC = points[2][0]
    DistC = points[2][2]*MagicDistanceFactor

    # Transform from Latitude/Longitude to ECEF coordinates.
    xA,yA,zA = pyproj.transform(lla, ecef, LonA, LatA, 0, radians=False)
    xB,yB,zB = pyproj.transform(lla, ecef, LonB, LatB, 0, radians=False)
    xC,yC,zC = pyproj.transform(lla, ecef, LonC, LatC, 0, radians=False)

    # Convert to numpy arrays.
    P1 = numpy.array([xA/1000.0, yA/1000.0, zA/1000.0])
    P2 = numpy.array([xB/1000.0, yB/1000.0, zB/1000.0])
    P3 = numpy.array([xC/1000.0, yC/1000.0, zC/1000.0])

    # Sphere intersection from Wikipedia + Stackoverflow.
    ex = (P2 - P1)/(numpy.linalg.norm(P2 - P1))
    i = numpy.dot(ex, P3 - P1)
    ey = (P3 - P1 - i*ex)/(numpy.linalg.norm(P3 - P1 - i*ex))
    ez = numpy.cross(ex,ey)
    d = numpy.linalg.norm(P2 - P1)
    j = numpy.dot(ey, P3 - P1)
    x = (pow(DistA,2) - pow(DistB,2) + pow(d,2))/(2*d)
    y = ((pow(DistA,2) - pow(DistC,2) + pow(i,2) + pow(j,2))/(2*j)) - ((i/j)*x)

    # Handle errors.
    if pow(DistA,2) - pow(x,2) - pow(y,2) < 0:
        return []
    z = numpy.sqrt(pow(DistA,2) - pow(x,2) - pow(y,2))

    lon,lat,altitude = pyproj.transform(ecef, lla, x*1000,y*1000,z*1000, radians=False)
    lon,lat,altitude = pyproj.transform(ecef, lla, x*1000,y*1000,z*1000, radians=True)

    triPt = P1 + x*ex + y*ey + z*ez
    lon,lat,altitude = pyproj.transform(ecef, lla, triPt[0]*1000,triPt[1]*1000,triPt[2]*1000, radians=False)

    return [lat,lon]


# This was copied from trilat_optproblem.py and changed for only three point
# combinations.
errors = []
results = []

# We run the algorithm for a subset of  3-combinations of inputs.
count = 0
for points in itertools.combinations(coords, 3):
    res = trilaterate(points)
    if len(res) == 0:
        continue
    results.append(res)
    errors.append(vincenty(real, [res[0],res[1]]))
    count += 1
    if count > 20: break

# Now we sort the output by the actual error and determine percentiles.
errors, results = zip(*sorted(zip(errors, results), key=operator.itemgetter(0), reverse=True))
median_entry = results[len(results)/2]
avg = reduce(lambda x, y: x + y, errors) / len(errors)
median = errors[len(errors)/2]
n90 = errors[int(len(errors)*0.1)]
print "%3d Median: %f Average: %f 90th-percentile: %f" % (3, median, avg, n90)
