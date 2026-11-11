# Engine wall condition bracket

Five legs of the same chamber are marched under different wall conditions on
the real nozzle contour, and the wall's effect on the gas is measured rather
than assumed. The legs are adiabatic, radiating, film-cooled, regenerative and
conjugate.

## What is shipped

- `record/` holds the distilled campaign. Each leg carries its per-field
  readings, its wall-duty series, its station breakdown and its steadiness
  residual, together with the digests of the campaign they came from. The
  record holds measurements, not conclusions.
- `gates.json` holds the judgement made from that record: each leg's heat
  duty, its settlement on the agreed screens, the chamber's shared mode and
  the tolerance of its period, and the cautions raised where a screen did not
  settle.
- `studies/` holds the three settlement studies: the calibration against a
  surrogate null, the detection-power study and the search-multiplicity study.
  Their numbers are quoted on the page and registered in the claim register.
- `../../docs/numerics/2d/engine-wall-bracket.md` is the page these
  artefacts publish into.

## Reproducing the evidence

Nothing here needs the marched fields. Every check below runs against the
shipped record.

Verify that the committed gate artefact still follows from the committed
record. The comparison is byte for byte when the two runs agree on the
encoding, and falls back to comparing the quantities themselves when a
platform's maths library rounds a borderline case the other way:

```bash
python3 cases/engine-wall-bracket/verify_committed_gate.py
```

Verify that the three shipped studies still follow from the committed record.
Every surrogate draw is keyed by a fixed seed, so the reproduction is exact:

```bash
python3 cases/engine-wall-bracket/verify_committed_studies.py
```

Republish both pages from the committed artefacts and confirm they do not
change:

```bash
python3 cases/engine-wall-bracket/publish_bracket_results.py \
  docs/numerics/2d/engine-wall-bracket.md
python3 cases/engine-wall-bracket/publish_bracket_results.py \
  docs/numerics/2d/conjugate-engine-walls.md
python3 cases/engine-wall-bracket/publish_settlement_tables.py \
  docs/numerics/2d/engine-wall-bracket.md \
  --power=cases/engine-wall-bracket/studies/power.json \
  --search=cases/engine-wall-bracket/studies/search.json \
  --calibration=cases/engine-wall-bracket/studies/calibration.json
git diff --stat -- docs/numerics/2d
```

## Marching the campaign again

`distill.py` rebuilds `record/` from a marched campaign, and it is the only
step that needs the fields. It writes the record and the digests that tie it
back to the campaign it came from:

```bash
python3 cases/engine-wall-bracket/distill.py --runs <campaign dir> --record record
```

Running a study against a marched campaign instead of the record is the same
command without the `--record` flag:

```bash
python3 cases/engine-wall-bracket/detection_power.py 800 <campaign dir> power.json
```

Drop `--record` to read the marched fields instead.

The judgement runs on the python standard library alone. It fits its
oscillation with householder reflections written out in the case itself, so
the verdict does not move with whichever maths library a runner happens to
ship, and the artefact's numbers are rounded before they are written.

## Gate surfaces

- `tools/case_gates.sh` runs this case in its `committed` mode, which judges
  the shipped record rather than marching a campaign.
- Continuous integration re-judges the record, re-runs the studies and
  republishes the pages, so a number on the page cannot drift away from the
  evidence beneath it.
