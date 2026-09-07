# Bonding descriptor matching: an intermittent CI failure and what is behind it

Investigation of the flaky `test_extractor` failure first seen on
[PR #46](https://github.com/gruenewald-lab/CGsmiles/pull/46)
([run 33726551486](https://github.com/gruenewald-lab/CGsmiles/actions/runs/33726551486),
2026-09-03). Written 2026-09-07.

Two branches are involved:

| branch | what it holds |
|---|---|
| `more_utils` | the extractor, the failing test, and the property tests added in `c33bb0a` |
| `fix/descriptor-matching-order` (off `master`) | the resolver fixes, `b007a7a` and `81037e9` |

The extractor does not exist on `master`, so the two halves of the story
necessarily live apart.

## 1. The symptom

`test_extractor[{[#C1]|6}.{#C1=[$]CCCC[$]}]` fails on roughly one CI run
in several with `AssertionError: Graphs not isomorphic`. The test resolves
a reference molecule, shuffles the atom keys and fragids with an unseeded
`random`, extracts fragments, writes them back out and re-resolves. Only
some shuffles fail, and the seed is not recorded, so a red run cannot be
replayed.

Reproduced locally at 2 of 60 seeds. The reference molecule is a linear
C24 alkane: 74 atoms in one piece. The round trip produced **78 atoms in
three pieces** — four extra hydrogens, which is two bonds that were never
made, with the dangling valences capped by `_adjust_hcount`.

Nothing raised. Nothing warned.

## 2. Root cause

Three separate defects, one per layer.

### 2.1 The extractor collapses distinct descriptor labels

`annotate_bonding_operators` ([`graph_utils.py:298`](cgsmiles/graph_utils.py))
gives every inter-fragment bond its own label — `>0`/`<0`, `>1`/`<1`, and
so on — which is unambiguous. `condense_fragments`
([`extractor.py:271`](cgsmiles/extractor.py)) then merges fragment
instances whose graphs are isomorphic, and the merge collapses those
labels. The output for the failing case was:

```
{[#C1B]([#C1B][#C1C][#C1D])[#C1A][#C1]}.{
  #C1  = ... [<0] ...
  #C1A = ... [<1] ... [>0] ...
  #C1B = ... [>1] ... [<1] ...
  #C1C = ... [<4] ... [>1] ...
  #C1D = ... [>4] ...
}
```

`#C1B` now carries both ends of label `1`, so three different meta edges
compete for it: `C1B#0–C1B#1`, `C1B#0–C1A`, and `C1B#1–C1C`.

The condensation criterion is graph isomorphism plus a guard against
neighbouring a branch point. The criterion it actually needs is that the
resulting fragment set still admits exactly one pairing up to
automorphism, which is a strictly stronger condition.

### 2.2 The resolver matched greedily, in insertion order

`edges_from_bonding_descrpt` walked `meta_graph.edges` in whatever order
the graph happened to store them, and `match_bonding_descriptors` returned
the first compatible pair it found. No backtracking, no notion of which
edge was more constrained. Tracing the failing string:

```
C1B#0 -- C1B#1   used ('>11', '<11')   available ['>11','<11'] / ['>11','<11']
C1B#0 -- C1A#4   DROPPED               available ['<11']       / ['<11','>01']
C1B#1 -- C1C#2   DROPPED               available ['>11']       / ['<41','>11']
C1C#2 -- C1D#3   used ('<41', '>41')
C1A#4 -- C1#5    used ('>01', '<01')
```

The `C1B#0–C1B#1` edge was the only one of the three with a *choice*, and
it went first, eating the two descriptors its neighbours needed. Note that
a correct pairing existed and was entirely forced — pure unit propagation
would have found it:

| demand | admissible pairs | count |
|---|---|---|
| `C1A#4 — C1#5` | `(>0, <0)` | 1 |
| `C1B#0 — C1A#4` | `(>1, <1)` | 1 |
| `C1B#1 — C1C#2` | `(<1, >1)` | 1 |
| `C1C#2 — C1D#3` | `(<4, >4)` | 1 |
| `C1B#0 — C1B#1` | `(>1, <1)`, `(<1, >1)` | 2 |

Assign the four forced ones, and the fifth collapses to one option too.

### 2.3 An unsatisfiable edge was dropped in silence

```python
except LookupError:
    continue
```

That is what turns a matching failure into a corrupted molecule instead of
an error message.

## 3. A fix that looks obvious and is wrong

The natural first move is to raise at the match site rather than
`continue`. **It breaks four to five existing tests**, and for a good
reason: the silent skip is load-bearing.

```
{[#SC3A]1[#SC3][#TP1]1}.{#SC3A=CCCC[!],#SC3=CCCC[$][!],#TP1=[$]O}
```

The meta graph is a triangle. `SC3A–SC3` consumes the two `[!]`,
`SC3–TP1` consumes the two `[$]`, and `TP1–SC3A` has nothing left. But the
edge is real: `#SC3` carries its `[$]` and its `[!]` on the *same* atom, so
once squashing merges that atom into `SC3A`, the merged atom already holds
the bond to `TP1`.

A meta edge can therefore be realised by a merge rather than by a bond, and
matching runs *before* `squash_atoms` ([`resolve.py:426-428`](cgsmiles/resolve.py)),
so at match time you cannot yet tell the two apart.

The consequence for design: the matching objective is to **maximise**
satisfied edges, not to satisfy all of them, and the "did every edge get
realised" question belongs in a check that runs *after* squashing.

## 4. What was implemented

### 4.1 `b007a7a` — serve the most constrained edge first

Build one demand per bond a meta edge asks for (an edge of order *n*
contributes *n*), compute each demand's admissible descriptor pairs, and
serve the demand with the fewest pairs first, recomputing after each
assignment. Only demands touching a bead that just gave up a descriptor
are recomputed. Ties keep the original edge order, so uncontended strings
resolve exactly as before.

Gated on `self.legacy`. Under `legacy=False`, `compatible()` truncates
descriptors to their first character ([`resolve.py:41-43`](cgsmiles/resolve.py)),
so every `<` matches every `>` regardless of label and there is nothing to
order by; that path keeps the original behaviour.

The descriptor scan was split out as `admissible_matches` — later
`_pair_options` — and `match_bonding_descriptors` reimplemented on top of
it, so there is one copy of the compatibility loop.

### 4.2 `81037e9` — backtrack over the pairs

Ordering is not sufficient. When several edges tie on how many pairs they
admit, the order cannot separate them, and the first pair tried for
whichever goes first may still be the wrong one:

```
{[#F0]([#F2])[#F1]}.{#F0=[<a]CC[>a],#F1=[>a]CCCC[>a],#F2=[<a]CCCC[>a]}
```

Both edges admit two pairs. Either ordering loses an edge; only trying the
other pair recovers it.

So when the greedy pass leaves an edge unmatched, the pairs are searched
with backtracking, most-constrained-first, deepening from zero skips up to
what greedy already achieved and keeping the greedy result if it cannot be
beaten. A step budget (`DESCRIPTOR_SEARCH_BUDGET = 20000`) bounds the
search and falls back to greedy on exhaustion.

The same commit replaced the greedy pass's linear scan for the next edge
with a lazy heap — see §5.2.

### 4.3 Regression tests

Six cases in `test_descriptor_matching_order`
([`test_molecule_resolve.py`](cgsmiles/tests/test_molecule_resolve.py)),
which separate all three versions:

| version | result |
|---|---|
| `master` | 4 of 6 fail |
| `b007a7a` (ordering) | 2 of 6 fail |
| `81037e9` (backtracking) | 6 pass |

## 5. Measurements

### 5.1 Correctness

3000 randomly generated shared-label molecules (3–7 beads, 2–3 bonding
sites, labels drawn from `{a}` or `{a, b}`). 1550 were genuinely
unsatisfiable by descriptor count and were discarded; the remaining 1450
were checked against a full backtracking oracle. Molecules left with at
least one unbonded meta edge:

| version | unbonded |
|---|---|
| `master` (insertion order) | 270 |
| `b007a7a` (ordering) | 119 |
| `81037e9` (backtracking) | **0** |

### 5.2 Performance

End-to-end resolution is unchanged — within noise, because parsing and
hydrogen rebuilding dominate. A 500-mer resolves in about 363 ms either
way. Timing `edges_from_bonding_descrpt` in isolation, in ms:

| case | `master` | ordering only | + backtracking | vs `master` |
|---|---|---|---|---|
| martini benzene | 0.024 | 0.050 | 0.040 | 1.6x |
| toluene (squash) | 0.027 | 0.051 | 0.046 | 1.7x |
| glucose (virtual edges) | 0.024 | 0.051 | 0.039 | 1.6x |
| PEO 50 (undirected) | 0.555 | 1.144 | 0.832 | 1.5x |
| PEO 200 (undirected) | 1.870 | 5.888 | 3.590 | 1.9x |
| PEO 200 (directed) | 2.000 | 5.893 | 3.102 | 1.6x |
| PEO 500 (undirected) | 3.944 | 24.297 | 9.741 | 2.5x |

Two things worth reading out of that table.

**The backtracking is not what costs.** It only runs when the greedy pass
fails, and the searches are small when it does; the last column is
*faster* than ordering-only nearly everywhere.

**The ordering commit shipped a quadratic bug.** PEO 500 went 3.9 → 24.3 ms
because the greedy pass picked its next edge with a linear scan over all
pending demands, giving O(E²) in beads. A lazy heap with stale-entry
rejection brought it back to 9.7 ms. Both versions now scale linearly:
0.83 / 3.59 / 9.74 for 50 / 200 / 500 beads. This only surfaced because
the step was timed on its own — it is invisible end to end.

The residual 1.6–2.5x constant factor is inherent. Knowing how many pairs
an edge admits means the scan can no longer stop at the first match, which
is what `match_bonding_descriptors` used to do.

## 6. What is still open

### 6.1 A post-resolution check (not implemented, belongs on `master`)

Since the eager raise is ruled out (§3), the check has to run after
`resolve()`: every meta edge of order > 0 must have at least one atomistic
bond crossing the two beads, **counting an atom that belongs to both** —
that is what credits a squash.

Prototyped and verified: zero false positives across the squash cases, the
crown ether, Martini glucose with virtual edges, benzene, the cyclohexane
bond-order linearisation and `|`-repeats; and on the broken string it names
exactly the two dropped edges, `(0, 4)` and `(1, 2)`. Dropping the
squash-credit term makes all three squash cases false-positive, so that
term is required.

This is the change that converts silent corruption into a precise
diagnostic, and it is independent of everything else.

### 6.2 Extractor self-verification (not implemented, belongs on `more_utils`)

Proving locally that a condensation stays unambiguous is the hard general
problem, but the extractor already holds the target molecule. After
condensing, resolve the string it just wrote, run the §6.1 check plus an
isomorphism check, and on failure back off — undo the last condensation for
the offending fragname, toward per-edge-unique labels — then retry.
Correct by construction, degrades to a verbose-but-right string, costs one
resolve per back-off round.

A narrower stopgap: never emit a definition carrying both `>k` and `<k` for
the same `k`. That is the shape in every observed failure, but it is a
heuristic — it does not cover three *different* fragments sharing one
label.

### 6.3 Ambiguity detection

Nothing above catches a string that resolves completely but to the wrong
isomer. Confirmed:

```
{[#TC5]1[#TC5][#TC5]1}.{#TC5=[$]cn[$]}
```

resolves with every meta edge realised and with a spurious C–C plus N–N
bond in a ring that should alternate. Catching it means enumerating
solutions past the first and comparing the *molecules* they produce —
"more than one solution" is the wrong test, since benzene
(`{#TC5=[$]cc[$]}`) has many pairings that all give benzene. Enumerate
under a budget, hash each realised molecule (element / charge / bond
order), and flag when two hashes differ.

This is also the machinery a canonical CGsmiles would need, since "reuse of
a fragment name is sound iff the pairing is unique up to automorphism" is
the same question.

### 6.4 Non-legacy mode

Untouched by design. Without labels the problem is one large equivalence
class and ambiguity is pervasive; `sample.py` also leans on non-legacy
label semantics for probabilities and would need review first.

### 6.5 Test flakiness

`c33bb0a` on `more_utils` adds `cgsmiles/tests/test_hypothesis.py`, which
generates the molecule, the fragment palette and the permutations with
hypothesis rather than an unseeded `random`, so failures replay and shrink.
Two properties hold today and guard the reader/writer; the two
extract/write/resolve round trips are marked `xfail(strict=False)` for the
§2.1 defect and reproduce it on roughly three runs in four.

Once §6.2 lands, those markers should come off and the tests should switch
back to `SETTINGS` so shrinking is enabled again.

The hand-written `test_extractor` should also move onto the generated
permutations, or at minimum seed `random`.

## 7. Reproducing the measurements

The correctness sweep and both benchmarks were run as throwaway scripts.
The shapes worth keeping:

- **Correctness** — generate random meta graphs, assign bodies with 2–3
  bonding sites and directed labels from a small alphabet, discard the ones
  a full backtracking oracle proves unsatisfiable, resolve the rest, and
  count meta edges with no crossing bond.
- **Isolated timing** — drive a `MoleculeResolver` up to just before
  `edges_from_bonding_descrpt` (set `meta_graph`, copy `atomname` to
  `fragname`, reset `molecule`, call `resolve_disconnected_molecule`), then
  time `deepcopy(template).edges_from_bonding_descrpt()` and subtract the
  cost of the `deepcopy` alone.
- **Version comparison** — `git show <ref>:cgsmiles/resolve.py` over the
  working file, run, restore.
