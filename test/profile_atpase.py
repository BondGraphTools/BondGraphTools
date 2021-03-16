from pyprof2calltree import convert, visualize
import logging
from cProfile import Profile
from atpase import atpase, bgt

FILENAME = "atpase_profile.kgrind"
LOG = "profile.log"
profiler = Profile()
logger = bgt.logger
logger.setLevel(logging.INFO)
handler = logging.FileHandler(filename=LOG, mode='w')
logger.addHandler(handler)

try:
    profiler.runcall(atpase)
except KeyboardInterrupt as ex:
    handler.flush()
    raise ex

convert(profiler.getstats(), FILENAME)
visualize(profiler.getstats())
