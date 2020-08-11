#! bin/python3

import pandas as pd
import numpy as np
import os,sys

FILE=sys.argv[1]

data = pd.read_csv(FILE, sep="\t", header=None)

positives = data.loc[data[4]=="+"].index.astype('int')
negatives = data.loc[data[4]=="-"].index.astype('int')

data.iloc[positives,2] = data.iloc[positives, 1]
data.iloc[negatives,1] = data.iloc[negatives, 2]

data.to_csv(str(FILE)+".TSS.tmp", sep="\t", index=False, header=False)
