import os
import re
import json
import shutil
import argparse
from tedana import workflows

parser = argparse.ArgumentParser()
parser.add_argument('--fmriprepDir', type=str, required=True)
parser.add_argument('--bidsDir', type=str, required=True)
parser.add_argument('--cores', type=int, default=1)
parser.add_argument('--subjectName', type=str, required=True)
args = parser.parse_args()

prep_data = args.fmriprepDir
bids_dir = args.bidsDir
subject = f"sub-{args.subjectName.upper()}"

# echos per sujeto
echo_files = []
for root, dirs, files in os.walk(prep_data):
    for f in files:
        if f.startswith(subject) and '_echo-' in f and f.endswith('_desc-preproc_bold.nii.gz'):
            echo_files.append(os.path.join(root, f))

echo_files = sorted(echo_files)

#  _echo- prefix
acq = re.search(r'(sub-[^/]+_.*)_echo-', os.path.basename(echo_files[0])).group(1)

json_files = []
for root, dirs, files in os.walk(bids_dir):
    for f in files:
        if (
            f.startswith(acq)
            and '_echo-' in f
            and f.endswith('_bold.json')
            and 'func' in root
            and 'derivatives' not in root
        ):
            json_files.append(os.path.join(root, f))

json_files = sorted(json_files)
print("JSON files usados:")
for jf in json_files:
    print("  ", jf)

echo_times = []
for jf in json_files:
    with open(jf) as h:
        echo_times.append(json.load(h)['EchoTime'])

# use a mask
brain_mask = os.path.join(
    os.path.dirname(echo_files[0]),
    f"{acq}_desc-brain_mask.nii.gz"
)

# Output dir
out_dir = os.path.join(os.path.abspath(os.path.dirname(prep_data)), "tedana", subject)

if os.path.isdir(out_dir):
    shutil.rmtree(out_dir)
os.makedirs(out_dir, exist_ok=True)

# Prefix non doble sub-sub-
prefix = f"{acq}_"

workflows.tedana_workflow(
    echo_files,
    echo_times,
    out_dir=out_dir,
    prefix=prefix,
    fittype="curvefit",
    tedpca="kic",
    verbose=True,
    gscontrol=None,
    n_threads=args.cores,
    mask=brain_mask
)
