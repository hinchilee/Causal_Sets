
import json
import sys
import os
import csv

path = ''
with open(path + 'max_distance.csv', 'a') as f:
    writer = csv.writer(f, lineterminator='\n')
    writer.writerow([1,2])
