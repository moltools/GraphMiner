#!/usr/bin/env python3
from .cli import cli
import argparse


def main():
    args = cli()
    print('Hi Giovi')
    exit(0)


if __name__ == "__main__":
    main()