"""Benchmark timing class"""
__author__ = "Keith Hughitt and Steven Christe"
__email__ = "keith.hughitt@nasa.gov"

import time
import datetime
import platform
import math
from optparse import OptionParser
from optparse import IndentedHelpFormatter

class BenchmarkTimer:
    """A simple benchmark timer class"""
    def __init__(self):
        """Instantiates a new BenchmarkTimer"""
        self.start_time = 0.0
        self.total_time = 0.0
        self.geom_time = 0.0
        self.test_num = 0
        
    def log(self, msg):
        """Print out a string with time taken for test"""
        end_time = time.time()
        
        self.test_num += 1
        
        elapsed_time = end_time - self.start_time
        
        self.total_time += elapsed_time
        self.geom_time += math.log(elapsed_time)
        
        output = '\t%d\t%f\t%s' % (self.test_num, elapsed_time, msg)
        print(output)

        self.start_time = time.time()

    def reset(self):
        """Reset the timer"""
        self.start_time = time.time()
        
    def parse_arguments(self):
        ''' Gets command-line arguments and handles validation '''
        parser = OptionParser("%prog [options]", 
                              formatter=IndentedHelpFormatter(4, 80))
        parser.add_option("-s", "--scale-factor", dest="scale_factor", 
                          type="int", help="factor to scale tests by", 
                          metavar="NUM", default=1)

        options, args = parser.parse_args() # pylint: disable=W0612

        return options
            
    def print_header(self, title):
        """Prints the output header containing system and time information"""
        header = ("|%s performance for Python %s (%s)\n"
                  "|\tOS_FAMILY=%s, OS=%s, ARCH=%s %s\n"
                  "|\t%s") % (
            title, platform.python_version(), platform.python_build()[0], 
            platform.system(), " ".join(platform.dist()), 
            platform.processor(), platform.machine(), 
            datetime.datetime.today().ctime()
        )

        #Display header information
        print(header)
            
    def print_summary(self):
        """Prints a summary of the test results"""
        geom_mean = math.exp(self.geom_time / self.test_num)
        summary = ("\t%f=Total Time,"
                   "\t%f=Geometric mean,"
                   "\t%d tests.") % (self.total_time, geom_mean, self.test_num)

        print(summary)

