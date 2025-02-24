# Copyright (c) 2025 Simula Research Laboratory
# SPDX-License-Identifier: GPL-3.0-or-later
from os import chdir, getcwd

import pytest


@pytest.fixture
def save_cwd():
    # Save the current working directory
    cwd = getcwd()

    # Yield control to the test function
    yield

    # Restore the original working directory after the test completes
    chdir(cwd)
