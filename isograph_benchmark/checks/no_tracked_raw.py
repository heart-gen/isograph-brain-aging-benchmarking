from __future__ import annotations

import subprocess
import sys


def main() -> None:
    proc = subprocess.run(
        ["git", "ls-files", "inputs/raw"],
        check=False,
        text=True,
        capture_output=True,
    )
    tracked = [line for line in proc.stdout.splitlines() if line.strip()]
    if tracked:
        print("Raw files are tracked but must remain ignored:", file=sys.stderr)
        for path in tracked:
            print(f"  {path}", file=sys.stderr)
        raise SystemExit(1)
    print("OK: no tracked files under inputs/raw")


if __name__ == "__main__":
    main()

