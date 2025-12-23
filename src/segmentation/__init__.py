"""
Segmentation module

Handles image segmentation tasks including SAM3-based segmentation
for roads, trees, and vegetation detection.
"""

from .sam3_segment import segment_image

__all__ = ["segment_image"]

