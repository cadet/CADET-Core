"""Compare the old and the fixed quadrature node counts of the DG geometry operators.

degree of the integrand  = axDispQuadDeg + geomFactorDegree + 2 * polyDeg
required number of nodes = ceil((degree + 1) / 2), since n Gauss nodes are exact up to 2n - 1
"""

import math

for name, geom in (("radial", 1), ("frustum", 2), ("smoothly varying", None)):
    print(f"--- {name} ---")
    print("polyDeg axDispQuadDeg   required   old   new")
    for p in range(1, 6):
        g = p if geom is None else geom
        for q in range(0, 6):
            degree = q + g + 2 * p
            required = math.ceil((degree + 1) / 2)
            old = (q + g + 2 * p + 1) // 2          # std::ceil of an integer division
            new = (q + g + 2 * p + 2) // 2
            flag = "   <-- was short" if old < required else ""
            if old != new or flag:
                print(f"{p:7d} {q:13d} {required:10d} {old:5d} {new:5d}{flag}")
    print()
