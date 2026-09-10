"""Check that simple rst tables in a file have separator rows at least as wide as their cells.

Usage: python check_rst_tables.py <file.rst>
"""

import io
import re
import sys

path = sys.argv[1]
lines = io.open(path, encoding='utf-8').read().splitlines()

for i, line in enumerate(lines):
    if not re.fullmatch(r'\s*=+(\s+=+)+\s*', line):
        continue
    # separator row; the content row is the next line (skip if this closes a table)
    if i + 1 >= len(lines):
        continue
    content = lines[i + 1]
    if re.fullmatch(r'\s*=+(\s+=+)+\s*', content):
        continue
    spans = [(m.start(), m.end()) for m in re.finditer(r'=+', line)]
    for j, (start, end) in enumerate(spans):
        cell = content[start:end] if len(content) > start else ''
        beyond = content[end:end + 1] if len(content) > end else ' '
        if beyond.strip():
            print(f"{path}:{i + 2}: column {j + 1} overflows its separator "
                  f"(width {end - start}, needs at least {len(content[start:].rstrip())})")
            print("    " + content)
            print("    " + line)
print("checked " + path)
