"""
PyPASA - Python implementation of PASA alignment assembly pipeline.

This package provides a modern Python implementation of the core PASA
alignment assembly functionality, focusing on speed and maintainability
while maintaining compatibility with the original Perl implementation.
"""

__version__ = "0.1.0"
__author__ = "PASA Development Team"

from pypasa.models import Alignment, Assembly, AlignmentSegment

__all__ = [
    "Alignment",
    "Assembly", 
    "AlignmentSegment",
]
