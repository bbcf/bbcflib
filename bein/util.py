# bein/util.py
# Copyright 2010, BBCF

# This file is part of bein.

# Bein is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your
# option) any later version.

# Bein is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.

# You should have received a copy of the GNU General Public License
# along with bein.  If not, see <http://www.gnu.org/licenses/>.

"""
:mod:`bein.util` -- Library of functions for bein
=================================================

.. module:: bein.util
   :platform: Unix
   :synopsis: Useful programs and utilities for bein
.. moduleauthor:: Fred Ross <madhadron at gmail dot com>

This module is a repository of assorted robust functions that have
been developed for bein in day to day use.  Much of it is focused on
analysis of high throughput sequencing data, but if you have useful
functions for a different domain, please contribute them.
"""

import pickle
import re
import sys
import os
import threading
from contextlib import contextmanager

from bein import *

# Basic utilities

# 'cat' is used only as an example; it is useless in Bein.
# @program
# def _cat(input_file):
#     return {'arguments': ['cat',input_file],
#             'return_value': None}

# def cat(ex, input_file, filename=None):
#     if filename == None:
#         filename = unique_filename_in()
#     _cat(ex, input_file, stdout=filename)
#     return filename

# def _cat_nonblocking(ex, input_file, filename=None):
#     if filename == None:
#         filename = unique_filename_in()
#     f = _cat.nonblocking(ex, input_file, stdout=filename)
#     class Future(object):
#         def __init__(self, f):
#             self.future = f
#         def wait(self):
#             self.future.wait()
#             return filename
#     return Future(f)

# cat.nonblocking = _cat_nonblocking


def pause():
    """Pause until the user hits Return."""
    print >> sys.stdout, "Paused. Hit a key to continue."
    sys.stdin.readline()
    return None

def first_n_lines(input_file, n, output_file = None):
    """Writes the first *n* lines of *input_file* to another file.

    If *output_file* is ``None``, then the output is written to a randomly named file.
    """
    if output_file == None:
        output_file = unique_filename_in()
    with open(input_file, 'r') as inf:
        with open(output_file, 'w') as outf:
            for i in xrange(n):
                l = inf.readline()
                outf.write(l)
    return output_file


@program
def touch(filename = None):
    """Equivalent to shell: ``touch filename``

    Returns *filename*.  If filename is ``None``, *filename* is set to
    a unique, random name.
    """
    if filename == None:
        filename = unique_filename_in()
    return {"arguments": ["touch",filename],
            "return_value": filename}

@program
def remove_lines_matching(pattern, filename):
    output_file = unique_filename_in()
    return {'arguments': ['awk',"""!/%s/ { print $0 > "%s" }""" % (pattern,output_file),
                          filename],
            'return_value': output_file}


@program
def md5sum(filename):
    """Calculate the MD5 sum of *filename* and return it as a string."""
    def parse_output(p):
        m = re.search(r'=\s*([a-f0-9A-F]+)\s*$',
                      ''.join(p.stdout))
        return m.groups()[-1] # in case of a weird line in LSF
    return {"arguments": ["openssl","md5",filename],
            "return_value": parse_output}



@program
def sleep(n):
    """Sleep for *n* seconds.  Returns *n*."""
    return {"arguments": ["sleep", str(n)],
            "return_value": n}


@program
def count_lines(filename):
    """Count the number of lines in *filename* (equivalent to ``wc -l``)."""
    def parse_output(p):
        m = re.search(r'^\s*(\d+)\s+' + filename + r'\s*$',
                      ''.join(p.stdout))
        return int(m.groups()[-1]) # in case of a weird line in LSF
    return {"arguments": ["wc","-l",filename],
            "return_value": parse_output}


@program
def split_file(filename, n_lines = 1000, prefix = None, suffix_length = 3):
    """Equivalent to Unix command ``split``.

    *filename* is the file to split.  Returns a list of the names of
    the new files created.  *n_lines* is the number of lines to put in
    each file, *prefix* is the file prefix to use (which is set to a
    unique, randomly chosen string if not specified), and
    *suffix_length* is the number of positions to use after the prefix
    to label the files.
    """
    if prefix == None:
        prefix = unique_filename_in()
    def extract_filenames(p):
        return [x for x in os.listdir('.') if x.startswith(prefix)]
    return {"arguments": ["split", "-a", str(suffix_length),
                          "-l", str(n_lines), filename, prefix],
            "return_value": extract_filenames}


def use_pickle(ex_or_lims, id_or_alias):
    """Loads *id_or_alias* as a pickle file and returns the pickled objects.

    *ex_or_lims* may be either an execution object or a MiniLIMS object.
    """

    if isinstance(ex_or_lims, MiniLIMS):
        lims = ex_or_lims
    elif isinstance(ex_or_lims, Execution):
        lims = ex_or_lims.lims
    else:
        raise ValueError("ex_or_lims must be a MiniLIMS or Execution.")

    f = lims.path_to_file(id_or_alias)
    with open(f) as q:
        d = pickle.load(q)
    return d

def background(fun, *args, **kwargs):
    """Run a function, but return a Future object instead of blocking.

    Instead of blocking, it starts the function in a separate thread,
    and returns an object which lets the user choose when to wait for
    the function by calling its wait() method.  wait() blocks its
    current thread until the function returns, then wait returns the
    value returned by the function.

        f = background(sqrt, 0)
        a = f.wait()

    is exactly equivalent to

        a = sqrt(0)

    except that in the first case, sqrt is run in a separate thread.

    The argument list after *fun* is exactly what you would pass to
    *fun* if you were calling it directly, including keyword
    arguments.
    """
    class Future(object):
        def __init__(self):
            self.return_value = None

        def wait(self):
            v.wait()
            return self.return_value
    future = Future()
    v = threading.Event()
    def g():
        future.return_value = fun(*args, **kwargs)
        v.set()
    a = threading.Thread(target=g)
    a.start()
    return(future)

def deepmap(f, st):
    """Map function *f* over a structure *st*.

    *f* should be a function of one argument.  *st* is a data
    structure consisting of lists, tuples, and dictionaries.  *f* is
    applied to every value in *st*.

    >>> deepmap(lambda x: x, 1)
    1
    >>> deepmap(lambda x: x, [1, 2, 3])
    [1, 2, 3]
    >>> deepmap(lambda x: x+1, {'a': [1,2], 'b': [3,4]})
    {'a': [2, 3], 'b': [4, 5]}
    >>> deepmap(lambda x: x, (1, 2, 3))
    (1, 2, 3)
    >>> deepmap(lambda x: x, {1: 2, 3: 4})
    {1: 2, 3: 4}
    >>> deepmap(lambda x: x, {1: (2, [3, 4, 5]), 2: {5: [1, 2], 6: (3, )}})
    {1: (2, [3, 4, 5]), 2: {5: [1, 2], 6: (3,)}}
    """
    if isinstance(st, list):
        return [deepmap(f, q) for q in st]
    elif isinstance(st, tuple):
        return tuple([deepmap(f, q) for q in list(st)])
    elif isinstance(st, dict):
        return dict([(k,deepmap(f, v)) for k,v in st.iteritems()])
    else:
        return f(st)

# Adding special file types


def add_pickle(execution, val, description="", alias=None):
    """Pickle *val*, and add it to the repository.

    add_pickle lets you dump almost any Python value to a file in the
    MiniLIMS repository.  It is useful to keep track of intermediate
    calculations.  *description* will be set as the pickle file's
    description.
    """
    if isinstance(description,dict): description = str(description)
    filename = unique_filename_in()
    with open(filename, 'wb') as f:
        pickle.dump(val, f)
    execution.add(filename, description=description, alias=alias)
    return filename

try:
    import pylab
    @contextmanager
    def add_figure(ex, figure_type='eps', description="", alias=None, figure_size=None):
        """Create a matplotlib figure and write it to the repository.

        Use this as a with statement, for instance::

            with add_figure(ex, 'eps') as fig:
                hist(a)
                xlabel('Random things I found')

        This will plot a histogram of a with the x axis label set, and
        write the plot to the repository as an EPS file.
        """
        if isinstance(description,dict): description = str(description)
        f = pylab.figure(figsize=figure_size)
        yield f
        filename = unique_filename_in() + '.' + figure_type
        f.savefig(filename)
        ex.add(filename, description=description, alias=alias)
except:
    print >>sys.stderr, "Could not import matplotlib.  Skipping add_figure."



try:
    import tables as h5
    @contextmanager
    def add_hdf5(ex, description='', alias=None):
        if isinstance(description,dict): description = str(description)
        h5filename = unique_filename_in()
        db = h5.openFile(h5filename, 'w', title=description)
        try:
            yield db
        finally:
            db.close()
            ex.add(db, description=description, alias=alias)
except:
    print >>sys.stderr, "PyTables not found.  Skipping."



if __name__ == "__main__":
    import doctest
    doctest.testmod()
