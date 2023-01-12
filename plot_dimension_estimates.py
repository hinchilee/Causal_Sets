import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('dimension_estimates.csv')
plt.errorbar(df['N'], df['mean'], yerr=df['sem'], fmt='o')
plt.xscale('log')
plt.xlabel('N')
plt.ylabel('Dimension')
plt.show()