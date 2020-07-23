#!/bin/env python
import re

f = open("matches")

for line in f:
    m = re.match(r"^(\d+)\s([^\s]+)\s([^$])$", line)
    if (m == None): continue
    print m.group(0)

f.close()
