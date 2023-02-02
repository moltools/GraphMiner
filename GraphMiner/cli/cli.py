#!/usr/bin/env python3
import argparse 
#import pandas as pd

def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, default=True, help="input path")
    return parser.parse_args()