"""Tests for funvip.src.validate_option preset resolution.

Guards the preset-handler fixes: the substring->tuple membership conversion, the
CLUSTER-CUTOFF -> cluster.cutoff fix, the ~-on-bool fixes, and the e-value
default/preset consolidation.
"""
import os

import pytest

import funvip
from funvip.src.validate_option import Option

PRESET_DIR = os.path.join(os.path.dirname(funvip.__file__), "preset")


def _load(preset_name):
    opt = Option()
    opt.update_from_preset(os.path.join(PRESET_DIR, f"{preset_name}.yaml"))
    return opt


def test_option_defaults():
    opt = Option()
    assert opt.cluster.cutoff == 0.95
    assert opt.cluster.evalue == 0.0001  # default e-value is 1e-4
    assert opt.solveflat is True
    assert opt.method.tcs is True


@pytest.mark.parametrize("preset", ["fast", "accurate"])
def test_preset_resolves_cleanly(preset):
    opt = _load(preset)
    # CLUSTER-CUTOFF now lands on cluster.cutoff (was silently discarded to evalue)
    assert opt.cluster.cutoff == 0.97
    # single clean e-value key resolves to 1e-4
    assert opt.cluster.evalue == 0.0001
    # ~-on-bool fixes: these resolve to real booleans, not -2
    assert opt.solveflat is True
    assert opt.cachedb is True
    # substring-collision fixes: these keep their proper defaults, not a collided value
    assert opt.allow_innertrimming is False
    assert opt.nosearchresult is False
    assert opt.method.search == "blast"
    assert opt.method.trim == "trimal"


def test_fast_and_accurate_differ_on_tree_and_confident():
    fast = _load("fast")
    accurate = _load("accurate")
    assert fast.method.tree == "fasttree"
    assert accurate.method.tree == "raxml"
    assert fast.confident is True
    assert accurate.confident is False  # was the 'fale' typo
