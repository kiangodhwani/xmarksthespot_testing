# Complete workstation handoff

This snapshot includes all three downloaded BigBrain volumes and all saved results from this workstation, alongside the source, tests and experiment notes. `ARTIFACT_MANIFEST.json` records SHA-256 checksums and sizes for every file under `data/` and `results/`.

Start with [APPROACH_PROGRESSION.md](APPROACH_PROGRESSION.md) for the experiment narrative and existing figures, and [OUTPUTS.md](OUTPUTS.md) for the registration-ready CSV release. The final result package is `results/exports/landmark_release_v1.zip`.

The snapshot preserves historical, duplicate, preliminary and incomplete artifacts too. In particular, `landmark_release_v1.incomplete` is an abandoned partial export, and the -30 degree boundary preparation was never matched. Use the final release manifest and summary to identify completed tests; do not count every directory as an independent experiment.

## Restore on another machine

Clone the repository, enter its directory, then create a fresh environment:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests
```

There were 34 passing tests at handoff. Resolved package versions are also recorded with the saved matching runs. A copied virtual environment is not portable; `.venv/`, Python caches and authentication credentials are deliberately excluded.

To verify the saved data and results:

```bash
python - <<'PY'
import hashlib, json
from pathlib import Path
manifest = json.loads(Path('ARTIFACT_MANIFEST.json').read_text())
for item in manifest['files']:
    path = Path(item['path'])
    assert path.stat().st_size == item['bytes'], path
    with path.open('rb') as stream:
        assert hashlib.file_digest(stream, 'sha256').hexdigest() == item['sha256'], path
print(f"Verified {len(manifest['files'])} artifacts")
PY
```

Saved provenance contains original workstation paths as historical records. Consult the experiment notes for relative paths and reproduction commands; do not assume those absolute paths exist after cloning. Generators generally refuse to overwrite saved results.

BigBrain data and derived imagery retain their attribution and CC BY-NC-SA terms; see [README.md](README.md). The source URLs and expected source-volume hashes are in `scripts/download_data.py`.
