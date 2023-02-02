#!/usr/bin/env python3

from sys import argv

def load_data():
    df = pd.read_csv(argv[1])
    print(df)
    return 'YEEHES'


