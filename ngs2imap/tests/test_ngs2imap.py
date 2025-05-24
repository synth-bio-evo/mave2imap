"""
Unit and regression test for the ngs2imap package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import ngs2imap


def test_ngs2imap_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "ngs2imap" in sys.modules
