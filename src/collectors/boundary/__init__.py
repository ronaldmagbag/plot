"""
Boundary detection module

Modular boundary collector with separate components for:
- PropertyLine: Actual property boundaries
- SetbackLine: Setback boundaries (to be defined)
- BuildableEnvelope: Buildable area (to be defined)
- Classifier: Classify property lines into front/rear/sides
- INSPIRE: Cadastral data handling
"""

from .models import PropertyLine, SetbackLine, BuildableEnvelope, PropertyLineSegment
from .collector import BoundaryCollector

__all__ = [
    "PropertyLine",
    "SetbackLine",
    "BuildableEnvelope",
    "PropertyLineSegment",
    "BoundaryCollector",
]

