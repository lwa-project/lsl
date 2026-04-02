"""
Standalone helper to upload old telemetry log files.

Invoked as a subprocess by the telemetry client:
    python -m lsl.misc._telemetry_send <cache_dir> <timeout>

Uses a file lock to ensure only one upload runs at a time across all
processes.
"""

import os
import sys
import json
import glob
import time
import fcntl
from urllib.request import urlopen
from urllib.parse import urlencode


def main(cache_dir, timeout):
    lock_path = os.path.join(cache_dir, 'send.lock')
    try:
        lock_fh = open(lock_path, 'w')
        fcntl.flock(lock_fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except (OSError, IOError):
        # Another process is already sending
        return

    try:
        filenames = glob.glob(os.path.join(cache_dir, 'telemetry_*.log'))
        if not filenames:
            return

        oldest = min(os.path.getmtime(f) for f in filenames)
        if time.time() - oldest < 7*86400:
            return

        payloads = []
        for filename in filenames:
            try:
                with open(filename, 'r') as lf:
                    payloads.append(json.loads(lf.read()))
            except Exception:
                pass

        if not payloads:
            return

        data = urlencode({'payload': json.dumps(payloads)})
        with urlopen('https://fornax.phys.unm.edu/telemetry/log.php',
                     data.encode(), timeout=timeout) as uh:
            status = uh.read()

        if status == b'':
            for filename in filenames:
                try:
                    os.unlink(filename)
                except OSError:
                    pass
    except Exception:
        pass
    finally:
        fcntl.flock(lock_fh, fcntl.LOCK_UN)
        lock_fh.close()


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]))
