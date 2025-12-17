#!/usr/bin/env python3


def test_py_import():
    try:
        import PymiRa

        return True, "Import test passed."
    except Exception as e:
        return False, f"Import of PymiRa failed due to {e}"
