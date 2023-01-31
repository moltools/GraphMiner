#!/usr/bin/env python3
from .cli import cli 


def main():
    args = cli()
    print("Hello Giovi!")
    exit(0)


if __name__ == "__main__":
    main()