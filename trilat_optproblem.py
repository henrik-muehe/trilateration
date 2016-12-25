import scipy
import scipy.optimize
from haversine import haversine
from vincenty import vincenty
import itertools
import random
from pprint import pprint
import operator

# Test data; these could be your inputs. This is actually scraped of a website
# and not generated. The website in question only supplies real positions for
# a subset of search results but we use those to validate accuracy.
real = [48.14696,11.73628]
coords = [11.580218296782743,48.13918807515703,11.477],[11.517138039698214,48.18768612414934,16.672],[11.5608050377351,48.153738136990604,12.889],[11.617346446831283,48.11674732163194,9.331],[11.699343887674134,48.049757497992914,10.99],[11.552912258276171,48.153474425856174,13.465],[11.53954333934099,48.07802969169499,16.294],[11.697149289797098,48.09768441519491,6.115],[11.552367644497043,48.090519951888496,14.844],[11.590195038225199,48.088195668989414,12.504],[11.508514050229243,48.1922851309257,17.418],[11.478152945885506,48.1040735748235,19.512],[11.702692394012457,48.20727395298781,7.053],[11.63383506324293,48.13446046433374,7.637],[11.518997884978273,48.22400790883851,18.022],[11.597246216512007,48.0816412935017,12.462],[11.663962187263063,48.07012805230099,9.952],[11.715425377826055,48.12508875839166,2.842],[11.694984339584895,48.05326530794745,10.705],[11.684659307026715,48.094214333519645,6.91],[11.606752133062525,48.166948860628295,9.746],[11.703320847275219,48.19492282616845,5.784],[11.655273933378785,48.12798187910728,6.295],[11.568127687109664,48.07484573419233,14.653],[11.517470445042838,48.156072308432634,16.076],[11.548151207567555,48.0824097832868,15.511],[11.519083704698566,48.158870724575216,15.979],[11.60935979193895,48.22311257449991,12.495],[11.528556605113065,48.07641081289575,17.092],[11.478680391889089,48.12322489032875,19.073],[11.604630589559145,48.10073011293635,10.906]

# This section allows manipulating the input data.
#
# Setting accuracy_in_meters to 10 drops accuracy to 10 meters, setting it to
# 100 drops accuracy to 100 meters etc..
# Setting random_meters to anything other than 0 adds a random number between
# -random_meters and +random_meters to each input coordinate.
accuracy_in_meters = 1
random_meters = 0
for c in coords:
    c[2] = int((c[2]*1000)/accuracy_in_meters)*accuracy_in_meters
    c[2] = c[2] + random.uniform(-random_meters, random_meters);

# Function to optimize. We cheaply stop the algorithm from converging outside of
# the allowed geo coordinates by making any such coordinates return a high error
# value.
def func(x, points):
    if x[0] > 90 or x[1] > 90 or x[0] < -90 or x[1] < -90:
        return 9999999999
    distsum = 0
    for c in points:
        dist = vincenty([x[0],x[1]], [c[1],c[0]]) - x[2]*c[2]
        distsum += dist*dist
    return distsum

# We run the algorithm for any number of points.
for c in range(3,len(coords)+1):
    errors = []
    results = []

    # We run the algorithm for any combination of c-many points as long as there
    # are no more than 20. Since we have ~30 input points, this will otherwise
    # take forever.
    count = 0
    for points in itertools.combinations(coords, c):
        res = scipy.optimize.fmin_powell(func, [0,0,0], args=(points,), xtol=0.0001, ftol=0.0001, disp=False)
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
    print "%3d Median: %f Average: %f 90th-percentile: %f" % (c, median, avg, n90)
