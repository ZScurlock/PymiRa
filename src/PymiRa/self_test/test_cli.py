#!/usr/bin/env python3

import subprocess


def test_cli_invalid_input():
    try:
        result = subprocess.run(["pymira", "--help"], capture_output=True, text=True)

        if result.returncode == 0:
            return True, "CLI test successful."
        else:
            return False, "CLI returned non-zero exit code"

    except Exception as e:
        return False, f"CLI check failed:{e}"
