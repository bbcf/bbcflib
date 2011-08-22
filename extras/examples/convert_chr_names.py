'''
This scrpt will take a genomic track file and convert
the chromosme names from the 'chr4' format to the 'chrIV'
format

Usage:

python convert_chr_names.py fileIN fileOUT
'''

import sys, re
from bbcflib import track

def integer_to_roman(input):
    ints = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
    nums = ('M',  'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I')
    if type(input) != type(1): raise TypeError, 'Expected integer, got "%s."' % type(input)
    if not 0 < input < 4000: raise ValueError, 'Argument must be between 1 and 3999.'
    result = ""
    for i in range(len(ints)):
       count = int(input / ints[i])
       result += nums[i] * count
       input -= ints[i] * count
    return result

def convert(chrom):
    match = re.search('([a-zA-Z]+)([0-9]+)', chrom)
    return match.group(1) + integer_to_roman(int(match.group(2)))

with track.load(sys.argv[1]) as old:
    with track.new(sys.argv[2]) as new:
        for chrom in old:
            new.write(convert(chrom), old.read(chrom))
