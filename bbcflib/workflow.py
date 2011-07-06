"""
========================
Module: bbcflib.workflow
========================

A wrapper for functions to extend them to full workflows for HTSStation.
"""

# Other modules #
from bein import *
from bein.util import background

################################################################################
class Workflow(object):
    def __init__(self, f):
        self.f = f
        self.__name__ = f.__name__
        self.__doc__ = f.__doc__

    def __call__(self, ex, group1, group2, options):
        return self.f(ex, group1, group2, options)

    def workflow(self, lims, job_id, config):
        # Change to scratch dir defined in config

        # Get Job from frontend
        job = None

        # Build options dictionary
        options = {}

        # Calculate all appropriate group pairings
        group_pairs = []

        def _ex(i,j):
            description = 'job%d-%d-%d' % (frontend_id, i, j)
            with execution(lims, description=description) as ex:
                self.f(ex, job.groups[i], job.groups[j], options)
            return ((i,j), ex.id)

        futures = [background(_ex, i, j) for i,j in group_pairs]
        execution_ids = dict([f.wait() for f in futures])

        # Send email report of run
        return execution_ids

#-----------------------------------#
# This code was written by the BBCF #
# http://bbcf.epfl.ch/              #
# webmaster.bbcf@epfl.ch            #
#-----------------------------------#
