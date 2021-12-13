#!/usr/bin/env python
import defopt, sys
from . import tools
from .tools import *
from inspect import getmembers, isfunction


def entry_point():
    defopt.run([k[1] for k in getmembers(tools,isfunction)], version=True)


if __name__ == "__main__":
    entry_point()