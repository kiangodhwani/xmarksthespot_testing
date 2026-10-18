# Continue on Linux

Open the cloned repository folder in PyCharm. Use Python 3.12 and create a new
virtual environment on Linux; the Mac `.venv` cannot be reused.

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m unittest discover -s tests
python scripts/download_data.py
```

Run these from the repository root, in order, on a fresh clone:

```bash
python scripts/extract_pair.py --y 8 --output results/pairs/coronal_y8_repeat
python scripts/transform_pair.py
python scripts/generate_random_pairs.py --seed 42 --count 5
python scripts/mask_pair.py
python scripts/compare_reference_mask.py
python scripts/try_hysteresis.py
python scripts/detect_candidates.py --output results/candidates/harris_y8_peakfix
python scripts/test_descriptors.py
python scripts/explain_descriptor.py
python scripts/run_matching_baseline.py
```

Some generators deliberately refuse to overwrite existing output directories.
Use the commands on a fresh results folder, or inspect existing outputs before
choosing a new output path and updating dependent commands. No GPU is required
for the current baseline. Images are saved as PNG files: open them in PyCharm.
Data, generated results and virtual environments are excluded from Git. The
three source volumes download from BigBrain and are checked against our recorded
SHA-256 values. Keep the Mac folder if you want the original manual annotations
and original generated artifacts; regeneration does not recover manual clicks.
BigBrain attribution/licensing remains as described in the original README.

## Current state and next step

Read AGENTS.md, PROJECT_PATH.md and FOUNDATIONS.md before making changes.
Preparation uses already aligned volumes, the y=8 mm coronal plane, and known
synthetic histology transforms. Alignment is approximate; we did not establish a
universal anatomical error tolerance. The 1% histology hysteresis mask is
provisional search support, not validated tissue segmentation.

Independent Harris detection, polar-gradient descriptions, mutual/ratio matching,
RANSAC and reference evaluation now run end to end. The histology-only control
retains 117 correct pairs within 2 mm. MRI/histology is weak: five of six cases
return no consensus; one retains four pairs within 2 mm. These are augmentations
of one slice, not independent anatomical test cases. The next task is to review
failures and test one stronger cross-modality representation/context strategy
using the same evaluator. Do not change the benchmark thresholds silently.

For Codex CLI, start it from this repository folder and say:
“Read AGENTS.md, LINUX_SETUP.md, PROJECT_PATH.md and FOUNDATIONS.md. Continue the
MRI–histology project from the completed matching baseline. Explain the next
experiment before implementing it, and keep outputs interpretable.”
