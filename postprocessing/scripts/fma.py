import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.append('../../tests_sofia')
from utils import *

df = pd.read_pickle('../data/test_data')
plot_fma_from_file(df['qx_tot'], df['qy_tot'], df['c1'], df['x0'], df['y0'], save_to='../plots/test.png')
