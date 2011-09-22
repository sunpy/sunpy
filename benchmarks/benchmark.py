"""Benchmark timing class"""
__author__ = "Keith Hughitt and Steven Christe"
__email__ = "keith.hughitt@nasa.gov"

import time
import datetime
import platform
import math
import contextlib
import inspect
from optparse import OptionParser
from optparse import IndentedHelpFormatter
from collections import defaultdict

def stdist(values):
    ln = float(len(values))
    avg = sum(values) / ln
    stdev = math.sqrt(1 / (ln - 1) * sum((v - avg) ** 2 for v in values))
    return avg, stdev


class StopWatch(object):
    def __init__(self):
        self._start = None
        self.elapsed = 0
    
    def start(self):
        self._start = time.time()
        return self
    
    def stop(self, *args):
        self.elapsed += time.time() - self._start
    
    __enter__ = start
    __exit__ = stop
    
    def _pausefun(self):
        self.stop()
        yield
        self.start()
    
    def paused(self):
        return contextlib.contextmanager(self._pausefun)()


class BenchmarkResult(object):
    def __init__(self, bench, results=None, scaled=None):
        if results is None:
            results = []
        self.results = results
        self.scaled = scaled
        self.bench = bench
    
    def append(self, other):
        self.results.append(other)

    @property
    def stdist(self):
        return stdist(self.results)
    
    @property
    def avg(self):
        return stdist(self.results)[0]
    
    @property
    def stdev(self):
        return stdist(self.results)[1]
    
    @property
    def total(self):
        return sum(self.results)
    
    
class BenchmarkSuite(object):
    def __init__(self, title, benchs=None):
        if benchs is None:
            benchs = []
        self.title = title
        self.benchs = benchs
        # print [v for k, v in inspect.getmembers(self) if k.startswith('bench')]
        self.benchs.extend(
            v for k, v in inspect.getmembers(self)
            if k.startswith('bench') and inspect.ismethod(v)
        )
        
        self.benchs.sort()
    
    @staticmethod
    def _run_bench(bench, times=1):
        bres = BenchmarkResult(bench)
        for _ in xrange(times):
            with StopWatch() as watch:
                bres.scaled = bench(watch)
            bres.append(watch.elapsed)
        return bres
    
    def run(self, times=1):
        results = [self._run_bench(bench, times) for bench in self.benchs]
        
        elapsed = sum(res.total for res in results)
        gelapsed = sum(math.log(res.total) for res in results)
        return results, elapsed, gelapsed
    
    def print_result(self, results, elapsed, gelapsed):
        print (
            "|%s performance for Python %s (%s)\n"
            "|\tOS_FAMILY=%s, OS=%s, ARCH=%s %s\n"
            "|\t%s" % (
                self.title, platform.python_version(), platform.python_build()[0], 
                platform.system(), " ".join(platform.dist()), 
                platform.processor(), platform.machine(), 
                datetime.datetime.today().ctime()
            )
        )
        
        for res in results:
            print '\t%d\t%f (%f)\t%s' % (
                self.benchs.index(res.bench) + 1, res.avg, res.stdev, res.scaled)
            
        gelapsed = math.exp(gelapsed / len(results))
        print (
            "\t%f=Total Time,"
            "\t%f=Geometric mean,"
            "\t%d tests." % (elapsed, gelapsed, len(results))
        )
    
    def main(self, argv):
        parser = OptionParser("%prog [options]", 
                              formatter=IndentedHelpFormatter(4, 80))
        parser.add_option("-s", "--scale-factor", dest="scale_factor", 
                          type="int", help="factor to scale tests by", 
                          metavar="NUM", default=1)
        parser.add_option("-t", "--times", dest="times", 
                          type="int", help="run each test NUM  times", 
                          metavar="NUM", default=1)

        options, args = parser.parse_args(argv) # pylint: disable=W0612
        
        self.scale = options.scale_factor        
        self.print_result(*self.run(options.times))
