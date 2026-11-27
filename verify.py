"""Exact checks for the weak-gate path criterion and the radialene formulas.
Run with Python 3.8+ using only the standard library.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from functools import reduce
from math import gcd, prod


def determinant(matrix):
    """Row-scaled integer Bareiss elimination, with row pivoting."""
    if not matrix:
        return Q(1)
    rows = []
    scales = []
    for row in matrix:
        scale = reduce(lambda a, b: a // gcd(a, b) * b, (v.denominator for v in row), 1)
        rows.append([int(v * scale) for v in row])
        scales.append(scale)
    n = len(rows)
    sign = 1
    previous = 1
    for k in range(n - 1):
        if rows[k][k] == 0:
            pivot = next((i for i in range(k + 1, n) if rows[i][k]), None)
            if pivot is None:
                return Q(0)
            rows[k], rows[pivot] = rows[pivot], rows[k]
            sign *= -1
        pivot_value = rows[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                rows[i][j] = (rows[i][j] * pivot_value - rows[i][k] * rows[k][j]) // previous
            rows[i][k] = 0
        previous = pivot_value
    return Q(sign * rows[-1][-1], prod(scales))


def cofactor(matrix, t, s):
    minor = [[v for j, v in enumerate(row) if j != s]
             for i, row in enumerate(matrix) if i != t]
    return (-1) ** (s+t) * determinant(minor)


def shifted(matrix, gates):
    out = [row.copy() for row in matrix]
    for i, value in gates.items():
        out[i][i] -= value
    return out


def green(matrix, s, t, gates):
    m = shifted(matrix, gates)
    return cofactor(m, t, s) / determinant(m)


def connected_subcubic_graphs(n):
    """Enumerate labelled edge masks and discard full permutation orbits."""
    edges = list(combinations(range(n), 2))
    edge_index = {e: i for i, e in enumerate(edges)}
    images = [[1 << edge_index[tuple(sorted((p[i], p[j])))] for i, j in edges]
              for p in permutations(range(n))]
    seen = set()
    for mask in range(1 << len(edges)):
        if mask in seen or not n - 1 <= bin(mask).count("1") <= 3*n//2:
            continue
        ids = [i for i in range(len(edges)) if mask & (1 << i)]
        neighbors = [set() for _ in range(n)]
        for i in ids:
            u, v = edges[i]
            neighbors[u].add(v)
            neighbors[v].add(u)
        if max(map(len, neighbors)) > 3:
            continue
        reached = {0}
        stack = [0]
        while stack:
            v = stack.pop()
            new = neighbors[v] - reached
            reached.update(new)
            stack.extend(new)
        if len(reached) == n:
            seen.update(sum(row[i] for i in ids) for row in images)
            yield [[int(j in nbrs) for j in range(n)] for nbrs in neighbors]


def chordless_weights(g, s, t):
    """Enumerate induced paths. Index each by the occupied internal sites."""
    weights = {}

    def extend(path, weight):
        if path[-1] == t:
            support = frozenset(path[1:-1])
            assert support not in weights
            weights[support] = weight
            return
        for u, val in enumerate(g[path[-1]]):
            if val and u not in path and all(g[u][v] == 0 for v in path[:-1]):
                extend(path + [u], weight * val)
    extend([s], Q(1))
    return weights


def cofactor_coefficients(matrix, s, t, sites):
    """Normalized cofactor values at 0/1 corners, then Boolean Möbius inversion."""
    det0 = determinant(matrix)
    coeff = []
    for mask in range(1 << len(sites)):
        gates = {site: 1 for i, site in enumerate(sites) if mask & (1 << i)}
        coeff.append(cofactor(shifted(matrix, gates), t, s) / det0)
    for i in range(len(sites)):
        for mask in range(len(coeff)):
            if mask & (1 << i):
                coeff[mask] -= coeff[mask ^ (1 << i)]
    return coeff


def submasks(mask):
    subset = mask
    while subset:
        yield subset
        subset = (subset - 1) & mask
    yield 0


def sign(x):
    return (x > 0) - (x < 0)


def classify_pair(matrix, g, s, t):
    sites = [i for i in range(len(g)) if i not in (s, t)]
    coeff = cofactor_coefficients(matrix, s, t, sites)
    paths = {sum(1 << sites.index(v) for v in support): w
             for support, w in chordless_weights(g, s, t).items()}
    minimal = {m: val for m, val in enumerate(coeff) if val and
               all(not coeff[k] for k in submasks(m) if k != m)}
    assert minimal == paths, (matrix, s, t, minimal, paths)
    counts = Counter()
    for config in product((0, 1, -1), repeat=len(sites)):
        active = sum(1 << i for i, val in enumerate(config) if val)
        signs = [prod(config[i] for i in range(len(sites)) if m & (1 << i))
                 for m in range(len(coeff))]
        predicted = {sign(w * signs[m]) for m, w in paths.items() if (m & active) == m}
        # On each active corner, read the first nonzero univariate coefficient.
        actual = set()
        for corner in submasks(active):
            along_ray = [0] * (len(sites) + 1)
            for m in submasks(corner):
                along_ray[bin(m).count("1")] += coeff[m] * signs[m]
            first = next((v for v in along_ray if v), 0)
            if first:
                actual.add(sign(first))
        assert predicted == actual, (matrix, s, t, config, predicted, actual)
        counts[len(actual)] += 1
    return counts


def census():
    print("n  connected-subcubic  nonsingular  zero-contact-pairs  protected  one-sign  mixed")
    total = Counter()
    for n in range(2, 7):
        graphs = nonsingular = pairs = 0
        counts = Counter()
        for matrix in connected_subcubic_graphs(n):
            graphs += 1
            det0 = determinant(matrix)
            if not det0:
                continue
            nonsingular += 1
            g = [[cofactor(matrix, t, s) / det0 for t in range(n)] for s in range(n)]
            for s, t in combinations(range(n), 2):
                if g[s][t] == 0:
                    pairs += 1
                    counts.update(classify_pair(matrix, g, s, t))
        print(n, graphs, nonsingular, pairs, counts[0], counts[1], counts[2], sep="  ")
        total.update(counts)
    print("Total site/sign configurations checked:", sum(total.values()))


def radialene(n, ring, spoke):
    matrix = [[0] * (2*n) for _ in range(2*n)]
    for i in range(n):
        j = (i+1) % n
        matrix[i][j] = matrix[j][i] = ring[i]
        matrix[i][n+i] = matrix[n+i][i] = spoke[i]
    return matrix


def arc_value(path, ring, spoke, shifts):
    n = len(ring)
    edges = []
    for i, j in zip(path, path[1:]):
        edge = i if (i+1) % n == j else j
        edges.append(Q(ring[edge], spoke[i] * spoke[j]))
    sites = path[1:-1]
    older = latest = 1
    for j in range(1, len(sites)):
        older, latest = latest, latest - edges[j]**2 * shifts[sites[j-1]] * shifts[sites[j]] * older
    return (-1)**(len(path)-1) * prod(edges) * prod(shifts[i] for i in sites) / latest


def equal_shift_checks():
    # Minimal polynomials of 1/(2 cos(pi/n)), constant term first.
    factors = {5: [-1, 1, 1], 7: [1, -1, -2, 1], 9: [1, 0, -3, -1]}
    checked = 0
    for n, factor in factors.items():
        matrix = radialene(n, [1]*n, [1]*n)
        for t in range(2, n//2 + 1):
            sites = [n+i for i in range(n) if i not in (0, t)]
            coeff = cofactor_coefficients(matrix, n, n+t, sites)
            polynomial = [Q(0)] * (len(sites)+1)
            for mask, value in enumerate(coeff):
                polynomial[bin(mask).count("1")] += value
            assert any(polynomial)
            for degree in range(len(polynomial)-1, len(factor)-2, -1):
                quotient = polynomial[degree] / factor[-1]
                offset = degree - len(factor) + 1
                for j, value in enumerate(factor):
                    polynomial[offset+j] -= quotient * value
            assert not any(polynomial)
            checked += 1
    print("Equal-shift algebraic factors from full cofactors (n=5,7,9):", checked)


def radialene_checks():
    checked = 0
    for n in range(4, 13):
        for nonuniform in (False, True):
            ring = [1 + ((2*i+1) % 3 if nonuniform else 0) for i in range(n)]
            spoke = [1 + ((i+2) % 4 if nonuniform else 0) for i in range(n)]
            matrix = radialene(n, ring, spoke)
            for t in range(2, n-1):
                paths = [list(range(t+1)), [0] + list(range(n-1, t-1, -1))]
                for signed in (False, True):
                    shifts = {i: Q((-1 if signed and i % 2 else 1) * (i+1), 91)
                              for i in range(n) if i not in (0, t)}
                    exact = green(matrix, n, n+t, {n+i: v for i, v in shifts.items()})
                    assert exact == sum(arc_value(p, ring, spoke, shifts) for p in paths)
                    checked += 1
    print("Full-matrix vs two-arc identities at rational parameter points:", checked)
    u = Q(1, 10)
    for q in range(2, 13, 2):
        n = q+3
        matrix = radialene(n, [1]*n, [1]*n)
        older = latest = 1
        for _ in range(2, q+1):
            older, latest = latest, latest - u*u*older
        x = u**q / latest
        gates = {n+i: u for i in range(3, n)}
        gates[n+1] = x
        val = green(matrix, n, n+2, gates)
        assert val == 0
        for offset in (-1, 1):
            gates[n+1] = x * (1+Q(offset, 10))
            assert green(matrix, n, n+2, gates) == offset*x/10
        print(f"[{n}]radialene  y={u}  x={x}  G_st={val}  +/-10% x: opposite signs")
    n = 5
    eq = green(radialene(n, [1]*n, [1]*n), n, n+2, {n+i: u for i in (1,3,4)})
    assert eq == Q(89,990)
    print(f"[5]radialene  x=y={u}  G_st={eq}")
    y = Q(9, 10)
    x = y*y/(1-y*y)
    assert green(radialene(5, [1]*5, [1]*5), 5, 7, {6: x, 8: y, 9: y}) == 0
    print(f"[5]radialene  y={y}  x={x}  G_st=0")


if __name__ == "__main__":
    census()
    radialene_checks()
    equal_shift_checks()
