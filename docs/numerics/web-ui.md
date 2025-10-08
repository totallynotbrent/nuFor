# Web UI skeleton

The first pass at a web presence is deliberately minimal: a tiny, dependency-free
HTTP server in the `nufor` CLI that publishes the last computed snapshot as JSON
so a real front end can be layered on top of it.

## Running

    nufor serve [port]        # default port 8060

The server binds to 127.0.0.1, runs a Sod shock tube to t = 0.2, and serves:

- `/` — a minimal HTML page that fetches the snapshot and plots density.
- `/api/result` — the snapshot as JSON.

## API

`GET /api/result` returns the conserved and derived fields on the cell centers:

    {
      "n": 300, "gamma": 1.4, "time": 0.2,
      "centers": [...],        // x coordinates of cell centers
      "rho": [...],            // density
      "m": [...],              // momentum
      "e": [...],              // specific total energy
      "u": [...],              // velocity (derived)
      "p": [...],              // pressure (derived)
    }

The page is intentionally plain; the plan is to hand this skeleton over to an
external design tool and have it build the polished interface against the same
data API. A live residual stream and the config editor are follow-ups.