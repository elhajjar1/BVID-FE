"""Timelines shared by the animation tests, built once per module."""

from __future__ import annotations

import pytest

from tests.viz.anim._anim_helpers import small_config


@pytest.fixture(scope="module")
def sweep_timeline():
    from bvidfe.viz.anim.scenes.energy_sweep import build_timeline

    return build_timeline(small_config(), n_energies=10)


@pytest.fixture(scope="module")
def stack_timeline():
    from bvidfe.viz.anim.scenes.damage_stack import build_timeline

    return build_timeline(small_config())
