#!/usr/bin/env python3
import argparse 


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, default=True, help="input path")
    return parser.parse_args()