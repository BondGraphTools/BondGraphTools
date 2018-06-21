from cProfile import Profile
from atpase import atpase

FILENAME = "./stats/atpase_profile.kgrind"
profiler = Profile()
profiler.runcall(atpase)

from pyprof2calltree import convert, visualize

convert(profiler.getstats(), FILENAME)
visualize(profiler.getstats())

