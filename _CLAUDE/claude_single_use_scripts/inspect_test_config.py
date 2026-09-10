"""Print the structure of a CADET test config json, compactly."""

import json
import sys


def walk(o, depth=0, maxdepth=3):
    if depth > maxdepth or not isinstance(o, dict):
        return
    for k, v in o.items():
        if isinstance(v, dict):
            tag = "dict"
        elif isinstance(v, list):
            tag = "list[" + str(len(v)) + "]"
        else:
            tag = repr(v)[:50]
        print("  " * depth + k + " : " + tag)
        walk(v, depth + 1, maxdepth)


walk(json.load(open(sys.argv[1])))
