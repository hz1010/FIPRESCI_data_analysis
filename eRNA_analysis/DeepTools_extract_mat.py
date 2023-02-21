#### main function and package are from DeepTools github (https://github.com/deeptools/deepTools)
import sys
import pandas as pd
import gzip
from collections import OrderedDict
import numpy as np
from copy import deepcopy
import pyBigWig
from deeptools import getScorePerBigWigBin
from deeptools import mapReduce
from deeptools.utilities import toString, toBytes, smartLabels
from deeptools.heatmapper_utilities import getProfileTicks
import argparse
import os
import multiprocessing
from deeptools.parserCommon import writableFile, numberOfProcessors
from deeptools._version import __version__
from deeptools import parserCommon
from deeptools import heatmapper
import deeptools.computeMatrixOperations as cmo
import deeptools.deepBlue as db
import matplotlib.pyplot as plt


con=["Mix","PolyT","Random"]
for i in range(len(con)):
	path=con[i]+"_unique_2k.mat.gz"
	hm_RNA = heatmapper.heatmapper()
	print(path)
	hm_RNA.read_matrix_file(path)
	mat_RNA=hm_RNA.matrix.matrix
	aa=np.array(mat_RNA)
	np.savetxt("mat_"+con[i]+".txt",aa)

