"""
SAM3 Segmentation Script for Plot Analysis

Segments roads (not in OSM), trees, and grasses from Mapbox satellite imagery.
This script runs in the cu126 environment.
"""

import sys
import os
import json
import argparse
from pathlib import Path
from typing import Dict, Any, List, Tuple, Optional
import numpy as np

# Add sam3 project root to path (same pattern as sam3_annotator.py and sam3_mask_generator.py)
# This allows importing sam3 modules even if sam3 is not installed as a package
sam3_root = Path("E:/3dhouse/sam3")
if sam3_root.exists():
    sys.path.insert(0, str(sam3_root))

# Add project root to path for imports
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

# Note: This script should be run in the cu126 environment where sam3 is installed
# The micromamba command activates the environment before running this script

try:
    # Import dependencies - these should be available in cu126
    from PIL import Image
    import torch
    
    # Import directly from sam3 submodules (same pattern as working scripts)
    from sam3.model_builder import build_sam3_image_model
    from sam3.model.sam3_image_processor import Sam3Processor
    
    # Debug info
    print(f"Python executable: {sys.executable}")
    print(f"PyTorch version: {torch.__version__}")
    print(f"CUDA available: {torch.cuda.is_available()}")
    
except ImportError as e:
    print(f"Error importing SAM3 dependencies: {e}")
    print(f"Python executable: {sys.executable}")
    print(f"Python path (first 5): {sys.path[:5]}")
    print("\nTroubleshooting:")
    print("1. Make sure you're running this script via: micromamba run -n cu126 python ...")
    print("2. Verify sam3 is accessible: micromamba run -n cu126 python -c 'from sam3.model_builder import build_sam3_image_model; print(\"OK\")'")
    print("3. Check if sam3 directory exists at E:/3dhouse/sam3")
    print("4. If not installed, install it: micromamba run -n cu126 pip install -e E:/3dhouse/sam3")
    import traceback
    traceback.print_exc()
    sys.exit(1)


def segment_image(
    image_path: str,
    output_path: str,
    bbox: Optional[Dict[str, float]] = None,
    center_lon: Optional[float] = None,
    center_lat: Optional[float] = None
) -> Dict[str, Any]:
    """
    Segment roads, trees, and grasses from satellite imagery using SAM3.
    
    Args:
        image_path: Path to input satellite image
        output_path: Path to save segmentation results JSON
        bbox: Optional bounding box in degrees {west, east, south, north} for lat/lon conversion
        center_lon: Optional center longitude for lat/lon conversion
        center_lat: Optional center latitude for lat/lon conversion
    
    Returns:
        Dictionary with segmented features (pixel coordinates if bbox not provided)
    """
    print(f"Loading image from {image_path}")
    image = Image.open(image_path).convert("RGB")
    width, height = image.size
    print(f"Image size: {width}x{height}")
    
    # Initialize SAM3 model
    print("Initializing SAM3 model...")
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using device: {device}")
    
    # Build model - use default HuggingFace model (load_from_HF=True)
    # This automatically handles BPE file location
    model = build_sam3_image_model(
        checkpoint_path=None,
        load_from_HF=True,
        device=device
    )
    
    # Create processor with confidence threshold
    processor = Sam3Processor(
        model,
        device=device,
        confidence_threshold=0.3
    )
    
    # Set image in processor (only once)
    print("Setting image in processor...")
    inference_state = processor.set_image(image)
    print("Image set successfully")
    
    # Segment different features
    results = {
        "roads": [],
        "trees": [],
        "grasses": [],
        "metadata": {
            "image_size": [width, height],
            "bbox": bbox,
            "center": [center_lon, center_lat] if center_lon is not None and center_lat is not None else None,
            "coordinate_system": "lat_lon" if bbox is not None else "pixel"
        }
    }
    
    # Store raw masks for visualization (for debugging)
    raw_masks_dict = {
        "roads": [],
        "trees": [],
        "grasses": []
    }
    
    # 1. Segment roads (not in OSM)
    print("Segmenting roads...")
    try:
        # Reset prompts and set image again for new class (following sam3_annotator pattern)
        processor.reset_all_prompts(inference_state)
        inference_state = processor.set_image(image)
        
        # Set text prompt and get output
        output = processor.set_text_prompt(
            state=inference_state,
            prompt="road"
        )
        
        masks = output.get("masks", torch.tensor([]))
        scores = output.get("scores", torch.tensor([]))
        
        print(f"Found {len(masks)} road masks")
        
        # Store ALL raw masks for visualization (regardless of confidence)
        for i in range(len(masks)):
            mask = masks[i].squeeze().cpu().numpy()
            raw_masks_dict["roads"].append(mask)
        
        # Only convert high-confidence masks to polygons
        for i in range(len(masks)):
            if isinstance(scores[i], torch.Tensor):
                score = scores[i].item()
            else:
                score = float(scores[i])
            
            if score > 0.3:  # Confidence threshold
                mask = masks[i].squeeze().cpu().numpy()
                # Convert mask to polygon coordinates
                if bbox is not None:
                    polygons = mask_to_polygons(mask, bbox, width, height)
                else:
                    polygons = mask_to_polygons_pixel(mask, width, height)
                for poly in polygons:
                    results["roads"].append({
                        "id": f"road_seg_{i}_{len(results['roads'])}",
                        "type": "segmented_road",
                        "geometry": {
                            "type": "Polygon",
                            "coordinates": [poly]
                        },
                        "confidence": float(score)
                    })
        print(f"Added {len(results['roads'])} road segments (confidence > 0.3) out of {len(masks)} total masks")
    except Exception as e:
        print(f"Error segmenting roads: {e}")
        import traceback
        traceback.print_exc()
    
    # 2. Segment trees
    print("Segmenting trees...")
    try:
        # Reset prompts and set image again for new class
        processor.reset_all_prompts(inference_state)
        inference_state = processor.set_image(image)
        
        output = processor.set_text_prompt(
            state=inference_state,
            prompt="tree"
        )
        
        masks = output.get("masks", torch.tensor([]))
        scores = output.get("scores", torch.tensor([]))
        
        print(f"Found {len(masks)} tree masks")
        
        # Store ALL raw masks for visualization (regardless of confidence)
        for i in range(len(masks)):
            mask = masks[i].squeeze().cpu().numpy()
            raw_masks_dict["trees"].append(mask)
        
        # Only convert high-confidence masks to polygons
        for i in range(len(masks)):
            if isinstance(scores[i], torch.Tensor):
                score = scores[i].item()
            else:
                score = float(scores[i])
            
            if score > 0.3:
                mask = masks[i].squeeze().cpu().numpy()
                if bbox is not None:
                    polygons = mask_to_polygons(mask, bbox, width, height)
                else:
                    polygons = mask_to_polygons_pixel(mask, width, height)
                for poly in polygons:
                    results["trees"].append({
                        "id": f"tree_seg_{i}_{len(results['trees'])}",
                        "geometry": {
                            "type": "Polygon",
                            "coordinates": [poly]
                        },
                        "confidence": float(score)
                    })
        print(f"Added {len(results['trees'])} tree segments (confidence > 0.3) out of {len(masks)} total masks")
    except Exception as e:
        print(f"Error segmenting trees: {e}")
        import traceback
        traceback.print_exc()
    
    # 3. Segment grasses/vegetation
    print("Segmenting grasses and vegetation...")
    try:
        # Reset prompts and set image again for new class
        processor.reset_all_prompts(inference_state)
        inference_state = processor.set_image(image)
        
        output = processor.set_text_prompt(
            state=inference_state,
            prompt="grass"
        )
        
        masks = output.get("masks", torch.tensor([]))
        scores = output.get("scores", torch.tensor([]))
        
        print(f"Found {len(masks)} grass masks")
        
        # Store ALL raw masks for visualization (regardless of confidence)
        for i in range(len(masks)):
            mask = masks[i].squeeze().cpu().numpy()
            raw_masks_dict["grasses"].append(mask)
        
        # Only convert high-confidence masks to polygons
        for i in range(len(masks)):
            if isinstance(scores[i], torch.Tensor):
                score = scores[i].item()
            else:
                score = float(scores[i])
            
            if score > 0.3:
                mask = masks[i].squeeze().cpu().numpy()
                if bbox is not None:
                    polygons = mask_to_polygons(mask, bbox, width, height)
                else:
                    polygons = mask_to_polygons_pixel(mask, width, height)
                for poly in polygons:
                    results["grasses"].append({
                        "id": f"grass_seg_{i}_{len(results['grasses'])}",
                        "geometry": {
                            "type": "Polygon",
                            "coordinates": [poly]
                        },
                        "confidence": float(score)
                    })
        print(f"Added {len(results['grasses'])} grass segments (confidence > 0.3) out of {len(masks)} total masks")
    except Exception as e:
        print(f"Error segmenting grasses: {e}")
        import traceback
        traceback.print_exc()
    
    # Save results
    output_path_obj = Path(output_path)
    
    # If output_path is a directory, create a default filename
    if output_path_obj.is_dir() or (not output_path_obj.suffix and not output_path_obj.exists()):
        # Create default filename based on input image name
        input_image_name = Path(image_path).stem
        output_path_obj = output_path_obj / f"{input_image_name}_sam3_results.json"
        output_path = str(output_path_obj)
    
    # Ensure parent directory exists
    output_path_obj.parent.mkdir(parents=True, exist_ok=True)
    
    print(f"Saving results to {output_path}")
    with open(output_path, "w") as f:
        json.dump(results, f, indent=2)
    
    # Save raw masks visualization (for debugging)
    masks_vis_path = output_path_obj.parent / f"{output_path_obj.stem}_masks.png"
    print(f"Creating raw masks visualization: {masks_vis_path}")
    create_raw_masks_visualization(
        raw_masks_dict=raw_masks_dict,
        output_path=str(masks_vis_path),
        width=width,
        height=height
    )
    
    # Create comparison visualization: raw masks vs polygons (for roads only)
    comparison_vis_path = output_path_obj.parent / f"{output_path_obj.stem}_roads_comparison.png"
    print(f"Creating roads comparison visualization: {comparison_vis_path}")
    create_roads_comparison_visualization(
        raw_road_masks=raw_masks_dict.get("roads", []),
        road_polygons=results.get("roads", []),
        output_path=str(comparison_vis_path),
        width=width,
        height=height,
        bbox=bbox
    )
    
    # Create and save visualization image (from polygons)
    vis_path = output_path_obj.with_suffix('.png')
    print(f"Creating polygon visualization image: {vis_path}")
    create_visualization_image(
        image=image,
        results=results,
        output_path=str(vis_path),
        width=width,
        height=height,
        bbox=bbox
    )
    
    return results


def create_raw_masks_visualization(
    raw_masks_dict: Dict[str, List[np.ndarray]],
    output_path: str,
    width: int,
    height: int
) -> None:
    """
    Create a visualization of raw SAM3 masks (before polygon conversion).
    This helps debug mask extraction issues.
    
    Args:
        raw_masks_dict: Dictionary with class names as keys and lists of mask arrays as values
        output_path: Path to save visualization image
        width: Image width
        height: Image height
    """
    # Define colors for each class (RGB)
    class_colors = {
        "roads": (128, 128, 128),      # Gray
        "trees": (0, 128, 0),          # Dark Green
        "grasses": (0, 255, 0),        # Bright Green
    }
    
    # Create colored mask image
    colored_mask = np.zeros((height, width, 3), dtype=np.uint8)
    
    # Draw each class with its color
    for class_name, color in class_colors.items():
        masks = raw_masks_dict.get(class_name, [])
        if not masks:
            continue
        
        print(f"  Drawing {len(masks)} raw {class_name} masks in color {color}")
        
        for idx, mask in enumerate(masks):
            try:
                # Ensure mask is 2D and boolean
                if mask.ndim > 2:
                    mask = mask.squeeze()
                
                # Convert to boolean if needed
                if mask.dtype != bool:
                    mask_bool = mask > 0.5
                else:
                    mask_bool = mask
                
                # Ensure mask matches image dimensions
                if mask_bool.shape != (height, width):
                    # Resize if needed
                    mask_img = Image.fromarray(mask_bool.astype(np.uint8) * 255)
                    mask_img = mask_img.resize((width, height), Image.NEAREST)
                    mask_bool = np.array(mask_img) > 127
                
                # Count pixels in mask
                pixel_count = mask_bool.sum()
                if pixel_count == 0:
                    print(f"    Warning: Mask {idx} has no pixels, skipping")
                    continue
                
                # Draw mask with class color (use OR operation to combine overlapping masks)
                colored_mask[mask_bool] = color
                print(f"    Drawn mask {idx}: {pixel_count} pixels")
            except Exception as e:
                print(f"    Error drawing mask {idx}: {e}")
                continue
    
    # Save visualization
    vis_image = Image.fromarray(colored_mask)
    vis_image.save(output_path, quality=95)
    print(f"  ✓ Saved raw masks visualization: {output_path}")


def create_roads_comparison_visualization(
    raw_road_masks: List[np.ndarray],
    road_polygons: List[Dict[str, Any]],
    output_path: str,
    width: int,
    height: int,
    bbox: Optional[Dict[str, float]] = None
) -> None:
    """
    Create a side-by-side comparison of raw road masks vs extracted polygons.
    This helps debug polygon extraction issues.
    
    Args:
        raw_road_masks: List of raw mask arrays from SAM3
        road_polygons: List of road polygon dictionaries from results
        output_path: Path to save comparison image
        width: Image width
        height: Image height
        bbox: Optional bounding box for coordinate conversion
    """
    if not raw_road_masks and not road_polygons:
        print("  No road masks or polygons to compare")
        return
    
    # Create a comparison image: left side = raw masks, right side = polygons
    comparison = np.zeros((height, width * 2, 3), dtype=np.uint8)
    
    # Left side: Draw raw masks (gray)
    road_color = (128, 128, 128)  # Gray
    print(f"  Drawing {len(raw_road_masks)} raw road masks on left side")
    for idx, mask in enumerate(raw_road_masks):
        try:
            # Ensure mask is 2D and boolean
            if mask.ndim > 2:
                mask = mask.squeeze()
            
            # Convert to boolean if needed
            if mask.dtype != bool:
                mask_bool = mask > 0.5
            else:
                mask_bool = mask
            
            # Ensure mask matches image dimensions
            if mask_bool.shape != (height, width):
                mask_img = Image.fromarray(mask_bool.astype(np.uint8) * 255)
                mask_img = mask_img.resize((width, height), Image.NEAREST)
                mask_bool = np.array(mask_img) > 127
            
            pixel_count = mask_bool.sum()
            if pixel_count == 0:
                print(f"    Warning: Raw mask {idx} has no pixels, skipping")
                continue
            
            # Draw on left side
            comparison[:height, :width][mask_bool] = road_color
            print(f"    Drawn raw mask {idx}: {pixel_count} pixels")
        except Exception as e:
            print(f"    Error drawing raw mask {idx}: {e}")
            continue
    
    # Right side: Draw polygons (gray)
    print(f"  Drawing {len(road_polygons)} road polygons on right side")
    for idx, road in enumerate(road_polygons):
        geometry = road.get("geometry", {})
        coordinates = geometry.get("coordinates", [])
        
        if not coordinates or len(coordinates) == 0:
            continue
        
        # Get polygon coordinates
        if len(coordinates) > 0:
            if isinstance(coordinates[0], list) and len(coordinates[0]) > 0:
                if isinstance(coordinates[0][0], list):
                    poly_coords = coordinates[0]
                else:
                    poly_coords = coordinates
            else:
                continue
        else:
            continue
        
        # Convert to numpy array for drawing
        try:
            points = np.array(poly_coords, dtype=np.float64)
            
            # Convert lat/lon to pixels if needed
            if bbox is not None:
                for i, point in enumerate(points):
                    lon, lat = point[0], point[1]
                    norm_x = (lon - bbox["west"]) / (bbox["east"] - bbox["west"])
                    norm_y = 1 - (lat - bbox["south"]) / (bbox["north"] - bbox["south"])
                    points[i] = [norm_x * width, norm_y * height]
            
            # Convert to int32 for drawing
            points = points.astype(np.int32)
            
            # Ensure points are within image bounds
            points[:, 0] = np.clip(points[:, 0], 0, width - 1)
            points[:, 1] = np.clip(points[:, 1], 0, height - 1)
            
            # Draw on right side (offset by width)
            try:
                import cv2
                colored_mask_right = np.zeros((height, width, 3), dtype=np.uint8)
                colored_mask_right = np.ascontiguousarray(colored_mask_right, dtype=np.uint8)
                points_list = [np.ascontiguousarray(points, dtype=np.int32)]
                cv2.fillPoly(colored_mask_right, points_list, road_color)
                comparison[:height, width:width*2] = np.maximum(
                    comparison[:height, width:width*2],
                    colored_mask_right
                )
                pixel_count = (colored_mask_right > 0).sum() // 3  # Divide by 3 for RGB
                print(f"    Drawn polygon {idx}: {pixel_count} pixels (cv2)")
            except ImportError:
                from PIL import ImageDraw
                temp_mask = Image.new('L', (width, height), 0)
                draw = ImageDraw.Draw(temp_mask)
                draw.polygon([tuple(p) for p in points], fill=255)
                temp_array = np.array(temp_mask)
                comparison[:height, width:width*2][temp_array > 0] = road_color
                pixel_count = (temp_array > 0).sum()
                print(f"    Drawn polygon {idx}: {pixel_count} pixels (PIL)")
        except Exception as e:
            print(f"    Warning: Could not draw polygon {idx}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    # Add labels
    try:
        from PIL import ImageDraw, ImageFont
        comparison_img = Image.fromarray(comparison)
        draw = ImageDraw.Draw(comparison_img)
        
        try:
            font = ImageFont.truetype("arial.ttf", 20)
        except:
            font = ImageFont.load_default()
        
        # Label left side
        draw.text((10, 10), "Raw Masks", fill=(255, 255, 255), font=font)
        draw.text((10, 35), f"{len(raw_road_masks)} masks", fill=(200, 200, 200), font=font)
        
        # Label right side
        draw.text((width + 10, 10), "Polygons", fill=(255, 255, 255), font=font)
        draw.text((width + 10, 35), f"{len(road_polygons)} polygons", fill=(200, 200, 200), font=font)
        
        # Draw divider line
        draw.line([(width, 0), (width, height)], fill=(255, 255, 255), width=2)
        
        comparison_img.save(output_path, quality=95)
    except Exception as e:
        # Fallback: save without labels
        comparison_img = Image.fromarray(comparison)
        comparison_img.save(output_path, quality=95)
    
    print(f"  ✓ Saved roads comparison: {output_path}")
    print(f"    Left: {len(raw_road_masks)} raw masks | Right: {len(road_polygons)} polygons")


def create_visualization_image(
    image: Image.Image,
    results: Dict[str, Any],
    output_path: str,
    width: int,
    height: int,
    bbox: Optional[Dict[str, float]] = None
) -> None:
    """
    Create a colored visualization image showing all segmented classes.
    
    Args:
        image: Original input image
        results: Segmentation results dictionary
        output_path: Path to save visualization image
        width: Image width
        height: Image height
        bbox: Optional bounding box for coordinate conversion
    """
    # Define colors for each class (RGB)
    class_colors = {
        "roads": (128, 128, 128),      # Gray
        "trees": (0, 128, 0),          # Dark Green
        "grasses": (0, 255, 0),        # Bright Green
    }
    
    # Create colored mask image (mask only, no original image)
    colored_mask = np.zeros((height, width, 3), dtype=np.uint8)
    
    # Draw each class with its color
    for class_name, color in class_colors.items():
        features = results.get(class_name, [])
        if not features:
            continue
        
        print(f"  Drawing {len(features)} {class_name} segments in color {color}")
        
        for feature in features:
            geometry = feature.get("geometry", {})
            coordinates = geometry.get("coordinates", [])
            
            if not coordinates or len(coordinates) == 0:
                continue
            
            # Get polygon coordinates (first ring of polygon)
            # GeoJSON Polygon format: coordinates is [[[lon, lat], ...]] or [[[x, y], ...]]
            # We need the first (and usually only) ring
            if len(coordinates) > 0:
                if isinstance(coordinates[0], list) and len(coordinates[0]) > 0:
                    if isinstance(coordinates[0][0], list):
                        # Polygon format: [[[x, y], ...]] - get the ring
                        poly_coords = coordinates[0]
                    else:
                        # Already flat list: [[x, y], ...]
                        poly_coords = coordinates
                else:
                    # Empty or invalid structure
                    continue
            else:
                continue
            
            # Convert to numpy array for drawing
            try:
                points = np.array(poly_coords, dtype=np.float64)
                
                # Check if coordinates are in pixel space or lat/lon
                if bbox is not None:
                    # Convert lat/lon back to pixel coordinates
                    for i, point in enumerate(points):
                        lon, lat = point[0], point[1]
                        # Normalize to [0, 1]
                        norm_x = (lon - bbox["west"]) / (bbox["east"] - bbox["west"])
                        norm_y = 1 - (lat - bbox["south"]) / (bbox["north"] - bbox["south"])
                        # Convert to pixel coordinates
                        points[i] = [norm_x * width, norm_y * height]
                
                # Convert to int32 for drawing
                points = points.astype(np.int32)
                
                # Ensure points are within image bounds
                points[:, 0] = np.clip(points[:, 0], 0, width - 1)
                points[:, 1] = np.clip(points[:, 1], 0, height - 1)
                
                # Draw filled polygon
                try:
                    import cv2
                    # Ensure colored_mask is contiguous and properly formatted
                    colored_mask = np.ascontiguousarray(colored_mask, dtype=np.uint8)
                    # Ensure points is a list of numpy arrays
                    points_list = [np.ascontiguousarray(points, dtype=np.int32)]
                    cv2.fillPoly(colored_mask, points_list, color)
                except ImportError:
                    # Fallback: use PIL to draw
                    from PIL import ImageDraw
                    # Create a temporary image for this polygon
                    temp_mask = Image.new('L', (width, height), 0)
                    draw = ImageDraw.Draw(temp_mask)
                    draw.polygon([tuple(p) for p in points], fill=255)
                    temp_array = np.array(temp_mask)
                    colored_mask[temp_array > 0] = color
                except Exception as cv_error:
                    # If OpenCV fails, fall back to PIL
                    print(f"    OpenCV error, using PIL fallback: {cv_error}")
                    from PIL import ImageDraw
                    temp_mask = Image.new('L', (width, height), 0)
                    draw = ImageDraw.Draw(temp_mask)
                    draw.polygon([tuple(p) for p in points], fill=255)
                    temp_array = np.array(temp_mask)
                    colored_mask[temp_array > 0] = color
                    
            except Exception as e:
                print(f"    Warning: Could not draw {class_name} segment: {e}")
                continue
    
    # Save mask-only visualization (no original image overlay)
    vis_image = Image.fromarray(colored_mask)
    vis_image.save(output_path, quality=95)
    print(f"  ✓ Saved visualization (mask only): {output_path}")


def mask_to_polygons_pixel(
    mask: np.ndarray,
    img_width: int,
    img_height: int
) -> List[List[List[float]]]:
    """
    Convert binary mask to polygon coordinates in pixel space.
    
    Args:
        mask: Binary mask array (H, W) - values should be 0-1 or boolean
        img_width: Image width in pixels
        img_height: Image height in pixels
    
    Returns:
        List of polygons, each as [[x, y], ...] in pixel coordinates
    """
    try:
        from skimage import measure
    except ImportError:
        try:
            import cv2
            mask_uint8 = (mask * 255).astype(np.uint8) if mask.max() <= 1.0 else mask.astype(np.uint8)
            contours, _ = cv2.findContours(mask_uint8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            
            polygons = []
            for contour in contours:
                if len(contour) < 3:
                    continue
                
                poly_coords = []
                for pt in contour:
                    x, y = float(pt[0][0]), float(pt[0][1])
                    poly_coords.append([x, y])
                
                if poly_coords[0] != poly_coords[-1]:
                    poly_coords.append(poly_coords[0])
                
                polygons.append(poly_coords)
            
            return polygons
        except ImportError:
            # Last resort: return bounding box
            y_indices, x_indices = np.where(mask > 0.5)
            if len(y_indices) == 0:
                return []
            
            min_x, max_x = float(x_indices.min()), float(x_indices.max())
            min_y, max_y = float(y_indices.min()), float(y_indices.max())
            
            return [[[min_x, min_y], [max_x, min_y], [max_x, max_y], [min_x, max_y], [min_x, min_y]]]
    
    try:
        if mask.max() > 1.0:
            mask = mask / 255.0
        
        # Ensure mask is numpy array with proper dtype
        mask = np.asarray(mask, dtype=np.float64)
        
        # find_contours finds contours at the 0.5 threshold level
        # This means it finds the boundary between True (>0.5) and False (<0.5) regions
        contours = measure.find_contours(mask, 0.5)
        polygons = []
        
        print(f"    Found {len(contours)} contours from mask (shape: {mask.shape}, dtype: {mask.dtype}, min: {mask.min():.3f}, max: {mask.max():.3f})")
        
        for idx, contour in enumerate(contours):
            if len(contour) < 3:
                print(f"      Contour {idx}: Skipped (too few points: {len(contour)})")
                continue
            
            # Convert contour to numpy array and ensure proper dtype
            contour = np.asarray(contour, dtype=np.float64)
            
            poly_coords = []
            for point in contour:
                # measure.find_contours returns (row, col) = (y, x) in image coordinates
                # row (y) is the vertical position, col (x) is the horizontal position
                row, col = float(point[0]), float(point[1])
                # Save as [x, y] for drawing (x is column/horizontal, y is row/vertical)
                poly_coords.append([col, row])
            
            if len(poly_coords) > 0 and poly_coords[0] != poly_coords[-1]:
                poly_coords.append(poly_coords[0])
            
            if len(poly_coords) >= 3:  # Need at least 3 points for a polygon
                polygons.append(poly_coords)
                print(f"      Contour {idx}: Added polygon with {len(poly_coords)} points")
            else:
                print(f"      Contour {idx}: Skipped (too few points after processing: {len(poly_coords)})")
        
        print(f"    Total polygons extracted: {len(polygons)}")
        return polygons
    except Exception as e:
        print(f"Warning: Error extracting polygons: {e}")
        import traceback
        traceback.print_exc()
        return []


def mask_to_polygons(
    mask: np.ndarray,
    bbox: Dict[str, float],
    img_width: int,
    img_height: int
) -> List[List[List[float]]]:
    """
    Convert binary mask to polygon coordinates in lat/lon.
    
    Args:
        mask: Binary mask array (H, W) - values should be 0-1 or boolean
        bbox: Bounding box in degrees
        img_width: Image width in pixels
        img_height: Image height in pixels
    
    Returns:
        List of polygons, each as [[lon, lat], ...]
    """
    try:
        from skimage import measure
    except ImportError:
        print("Warning: scikit-image not available, using simplified polygon extraction")
        # Fallback: use OpenCV if available, or return bounding box
        try:
            import cv2
            # Convert mask to uint8
            mask_uint8 = (mask * 255).astype(np.uint8) if mask.max() <= 1.0 else mask.astype(np.uint8)
            contours, _ = cv2.findContours(mask_uint8, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
            
            polygons = []
            for contour in contours:
                if len(contour) < 3:
                    continue
                
                poly_coords = []
                for pt in contour:
                    x, y = pt[0][0], pt[0][1]
                    # Normalize to [0, 1]
                    norm_x = x / img_width
                    norm_y = y / img_height
                    
                    # Convert to lat/lon
                    lon = bbox["west"] + norm_x * (bbox["east"] - bbox["west"])
                    lat = bbox["south"] + (1 - norm_y) * (bbox["north"] - bbox["south"])
                    
                    poly_coords.append([lon, lat])
                
                # Close polygon
                if poly_coords[0] != poly_coords[-1]:
                    poly_coords.append(poly_coords[0])
                
                polygons.append(poly_coords)
            
            return polygons
        except ImportError:
            # Last resort: return bounding box of mask
            y_indices, x_indices = np.where(mask > 0.5)
            if len(y_indices) == 0:
                return []
            
            min_x, max_x = x_indices.min(), x_indices.max()
            min_y, max_y = y_indices.min(), y_indices.max()
            
            # Create bounding box polygon
            norm_x1, norm_x2 = min_x / img_width, max_x / img_width
            norm_y1, norm_y2 = min_y / img_height, max_y / img_height
            
            lon1 = bbox["west"] + norm_x1 * (bbox["east"] - bbox["west"])
            lon2 = bbox["west"] + norm_x2 * (bbox["east"] - bbox["west"])
            lat1 = bbox["south"] + (1 - norm_y1) * (bbox["north"] - bbox["south"])
            lat2 = bbox["south"] + (1 - norm_y2) * (bbox["north"] - bbox["south"])
            
            return [[[lon1, lat1], [lon2, lat1], [lon2, lat2], [lon1, lat2], [lon1, lat1]]]
    
    try:
        # Normalize mask to 0-1 range if needed
        if mask.max() > 1.0:
            mask = mask / 255.0
        
        # Ensure mask is numpy array with proper dtype
        mask = np.asarray(mask, dtype=np.float64)
        
        # Find contours
        contours = measure.find_contours(mask, 0.5)
        
        polygons = []
        for contour in contours:
            if len(contour) < 3:
                continue
            
            # Convert contour to numpy array and ensure proper dtype
            contour = np.asarray(contour, dtype=np.float64)
            
            # Convert pixel coordinates to lat/lon
            poly_coords = []
            for point in contour:
                # measure.find_contours returns (row, col) = (y, x) in image coordinates
                # row (y) is the vertical position, col (x) is the horizontal position
                row, col = float(point[0]), float(point[1])
                # Normalize to [0, 1]
                # col (x) maps to longitude (horizontal)
                # row (y) maps to latitude (vertical, but inverted - top is north)
                norm_x = col / img_width
                norm_y = row / img_height
                
                # Convert to lat/lon
                # x (col) -> longitude (west to east)
                lon = bbox["west"] + norm_x * (bbox["east"] - bbox["west"])
                # y (row) -> latitude (south to north, but image y=0 is top, so invert)
                lat = bbox["south"] + (1 - norm_y) * (bbox["north"] - bbox["south"])
                
                poly_coords.append([lon, lat])
            
            # Close polygon
            if len(poly_coords) > 0 and poly_coords[0] != poly_coords[-1]:
                poly_coords.append(poly_coords[0])
            
            if len(poly_coords) >= 3:  # Need at least 3 points for a polygon
                polygons.append(poly_coords)
        
        return polygons
    except Exception as e:
        print(f"Error converting mask to polygons: {e}")
        import traceback
        traceback.print_exc()
        return []


def main():
    parser = argparse.ArgumentParser(description="Segment roads, trees, and grasses from satellite imagery")
    parser.add_argument("--image", required=True, help="Path to input satellite image")
    parser.add_argument("--output", required=True, help="Path to output JSON file")
    parser.add_argument("--bbox", default=None, help="Optional bounding box as JSON: '{\"west\":...,\"east\":...,\"south\":...,\"north\":...}' for lat/lon conversion")
    parser.add_argument("--center-lon", type=float, default=None, help="Optional center longitude for lat/lon conversion")
    parser.add_argument("--center-lat", type=float, default=None, help="Optional center latitude for lat/lon conversion")
    
    args = parser.parse_args()
    
    # Parse bbox JSON if provided
    bbox = None
    if args.bbox:
        try:
            bbox = json.loads(args.bbox)
        except json.JSONDecodeError:
            print(f"Error: Invalid bbox JSON: {args.bbox}")
            sys.exit(1)
    
    # Run segmentation
    try:
        results = segment_image(
            args.image,
            args.output,
            bbox,
            args.center_lon,
            args.center_lat
        )
        print(f"Segmentation complete. Found:")
        print(f"  - {len(results['roads'])} road segments")
        print(f"  - {len(results['trees'])} tree segments")
        print(f"  - {len(results['grasses'])} grass/vegetation segments")
    except Exception as e:
        print(f"Error during segmentation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

