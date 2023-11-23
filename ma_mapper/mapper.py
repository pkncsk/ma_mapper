import pysam
import math
import numpy as np
import pandas as pd
#-> BAM file overlay 
def normal_func(x):
        return math.exp(-x * x / (2 * sigma2))

def normal_array(width=1, sigma=1, odd=0):
    ''' Returns an array of the normal distribution of the specified width '''
    sigma2 = float(sigma) ** 2
    # width is the half of the distribution
    values = list(map(normal_func, range(-width, width + odd)))
    values = np.array(values)

    return values

def bam_mapper(bam_file, metadata, offset = 5, probe_length = 100, smoothing_length = 100):
	bam_file = pysam.AlignmentFile(signal_file, "rb")
	