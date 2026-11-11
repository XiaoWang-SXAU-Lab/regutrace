# HTTP API plan (Phase 3)

ReguTrace starts as a CLI + SQLite backend. The server will expose the same
evidence store to the rsWGCNA dashboard and other clients.

Planned endpoints (versioned `/api/v1`):

```text
GET  /api/v1/health
GET  /api/v1/evidence?tf=&target=&level=&type=&source=
GET  /api/v1/tf/{id}/targets
GET  /api/v1/target/{id}/regulators
POST /api/v1/scan
POST /api/v1/import
GET  /api/v1/data-packs
```

Design rules:

- Return evidence records with level/provenance; never collapse A/B/C/D into a
  single opaque score.
- Keep predictions (D) separate from experimental evidence (A/B/C).
- Bind to localhost by default; the backend is local-first.
- The server and CLI share the same SQLite schema and importer.

The API will be implemented after Phase 1 (CLI + SQLite) and Phase 2 (imports +
benchmarks) are stable.
