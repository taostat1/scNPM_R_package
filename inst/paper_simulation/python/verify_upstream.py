"""Check AJGM source files against the source used for the manuscript."""
import hashlib
import json
from pathlib import Path
import sys

root = Path(sys.argv[1])
expected = json.loads(Path(__file__).with_name('upstream_AJGM_sha256.json').read_text())
for name, digest in expected.items():
    actual = hashlib.sha256((root / name).read_bytes().replace(b'\r\n', b'\n')).hexdigest()
    if actual != digest:
        raise SystemExit('Source differs: ' + name)
print('AJGM source fingerprints match.')
