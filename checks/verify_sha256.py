"""Check every file of this repository against SHA256SUMS.txt, or rewrite SHA256SUMS.txt and MANIFEST.md after a
deliberate change (any platform, Python 3, standard library only).

usage: python3 checks/verify_sha256.py           check (exit status 0 when every listed file is present and matches)
       python3 checks/verify_sha256.py --write   rewrite SHA256SUMS.txt and MANIFEST.md from the files as they are now
                                                 (for example after adding the DOI to README.md and CITATION.cff)
The script finds the repository's root itself. It skips .git/ and __pycache__/ and never lists MANIFEST.md or
SHA256SUMS.txt; when checking, files that SHA256SUMS.txt does not list (a regenerated solution.npz) are ignored."""
import hashlib
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKIP_DIRS = ('.git', '__pycache__')
SKIP_FILES = ('MANIFEST.md', 'SHA256SUMS.txt')


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def check():
    bad = ok = 0
    with open(os.path.join(ROOT, 'SHA256SUMS.txt'), encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            want, rel = line[:64], line[66:]
            path = os.path.join(ROOT, *rel.split('/'))
            if not os.path.isfile(path):
                print('MISSING  %s' % rel)
                bad += 1
            elif digest(path) != want:
                print('CHANGED  %s' % rel)
                bad += 1
            else:
                ok += 1
    print('%d files match, %d missing or changed' % (ok, bad))
    return 1 if bad else 0


def write():
    files = []
    for d, ds, fs in os.walk(ROOT):
        ds[:] = sorted(x for x in ds if x not in SKIP_DIRS)
        for f in sorted(fs):
            rel = os.path.relpath(os.path.join(d, f), ROOT).replace(os.sep, '/')
            if rel not in SKIP_FILES:
                files.append(rel)
    files.sort(key=lambda r: (r.count('/') > 0, r))
    rows = [(r, os.path.getsize(os.path.join(ROOT, r)), digest(os.path.join(ROOT, r))) for r in files]
    with open(os.path.join(ROOT, 'SHA256SUMS.txt'), 'w', encoding='utf-8', newline='\n') as fh:
        for r, n, h in rows:
            fh.write('%s  %s\n' % (h, r))
    with open(os.path.join(ROOT, 'MANIFEST.md'), 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('# Manifest\n\nEvery file of this repository, with its size in bytes and its SHA-256 (the same digests as '
                 '`SHA256SUMS.txt`, which `python3 checks/verify_sha256.py` or `sha256sum -c SHA256SUMS.txt` checks from '
                 'the repository\'s root). `MANIFEST.md` and `SHA256SUMS.txt` are not listed.\n\n')
        fh.write('| File | Bytes | SHA-256 |\n|---|---:|---|\n')
        for r, n, h in rows:
            fh.write('| `%s` | %d | `%s` |\n' % (r, n, h))
    print('wrote SHA256SUMS.txt and MANIFEST.md: %d files, %d bytes' % (len(rows), sum(n for _, n, _ in rows)))
    return 0


if __name__ == '__main__':
    sys.exit(write() if sys.argv[1:] == ['--write'] else check())
