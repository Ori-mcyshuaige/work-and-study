import os, sys, json, re, time, datetime, pickle
import matplotlib as mpl
import numpy as np, pandas as pd, matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as scs

from copy import copy, deepcopy
from glob import glob
from scipy.spatial.distance import pdist, squareform
from collections import defaultdict, OrderedDict
from functools import reduce

from .Variable import *

import warnings
warnings.filterwarnings('ignore')

#mpl.rcParams['font.sans-serif'] = ['simhei']
#np.random.seed(0)

def save_writer(writer):
    writer.save(); writer.close()