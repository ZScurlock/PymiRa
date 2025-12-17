#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  5 17:24:58 2025

@author: zac
"""
from PymiRa.self_test import test_import, test_bwt_align, test_cli


def run_self_test():
    results = []

    tests = [
        ("Importing package", test_import.test_py_import),
        ("CLI operations", test_cli.test_cli_invalid_input),
        ("Alignment", test_bwt_align.test_end_to_end),
    ]

    for name, test_func in tests:
        try:
            res, msg = test_func()
        except Exception as e:
            res, msg = False, f"Exception during test:{e}"
        results.append((name, res, msg))

    all_passed = all(res for _, res, _ in results)

    return all_passed, results


run_self_test()
