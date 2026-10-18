"""Download BigBrain inputs, checking recorded SHA-256 values."""
from pathlib import Path
import hashlib
import os
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
BASE = 'https://ftp.bigbrainproject.org/bigbrain-ftp/BigBrainRelease.2015/'
FILES = [
    ('mri_to_blk_aligned.nii.gz', BASE+'3D_MRI/nii/mri_to_blk_aligned.nii.gz', '098edc7b3a787890c39f12fb8fe31e6eaca703b341fd60ef6c653e0e269cf691'),
    ('full8_200um_optbal.nii.gz', BASE+'3D_Volumes/Histological_Space/nii/full8_200um_optbal.nii.gz', '69f19eb044d23a09c0c316e8ce985fcd6822914d19ac0ec36e72aa8e03026bef'),
    ('full_cls_200um_mask.nii.gz', 'https://bigbrain-ftp.loris.ca/bigbrain-ftp/BigBrainRelease.2015/3D_Classified_Volumes/Histological_Space/nii/full_cls_200um_mask.nii.gz', 'a5f55a76a5a01be256c4283a9f89b01e08afda8c3bfabc61e066b5118b9b24a7'),
]

def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def main():
    destination = ROOT/'data/bigbrain'
    destination.mkdir(parents=True, exist_ok=True)
    for name, url, expected in FILES:
        path = destination/name
        if path.exists():
            if digest(path) != expected:
                raise ValueError(f'Existing checksum mismatch: {path}; inspect before replacing')
            print(f'Verified existing {name}', flush=True)
            continue
        temporary = destination/(name+'.partial')
        print(f'Downloading {name}', flush=True)
        with urllib.request.urlopen(url, timeout=120) as source, temporary.open('wb') as output:
            while chunk := source.read(1024*1024):
                output.write(chunk)
        if digest(temporary) != expected:
            raise ValueError(f'Download checksum mismatch: {temporary}; destination unchanged')
        os.replace(temporary, path)
        print(f'Verified {name}', flush=True)

if __name__ == '__main__':
    main()
