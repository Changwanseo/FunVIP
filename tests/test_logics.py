"""Tests for funvip.src.logics — colour and NaN helpers."""
import numpy as np

from funvip.src.logics import isnan, isvalidcolor


def test_isvalidcolor_is_case_insensitive():
    # Regression: capitalized CSS names crashed option validation before the fix.
    assert isvalidcolor("Red")
    assert isvalidcolor("red")
    assert isvalidcolor("DarkGreen")


def test_isvalidcolor_hex():
    assert isvalidcolor("#aabbcc")
    assert isvalidcolor("#AABBCC")
    assert isvalidcolor("#abc")


def test_isvalidcolor_rejects_garbage():
    assert not isvalidcolor("notacolor")
    assert not isvalidcolor("#xyz123")


def test_isnan_no_longer_nameerrors():
    # Regression: isnan referenced np without importing numpy.
    assert isnan(float("nan")) is True
    assert isnan(np.float64("nan")) is True
    assert isnan(1.0) is False
    assert isnan("not a number") is False
