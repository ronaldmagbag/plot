
53.039271608492534, -1.2006180639806219


- Generation Plot Tests
python cli.py generate --lat 50.82991 --lon -0.266998 --output output\plot_adur1.json
python cli.py generate --lat 50.82981237 --lon -0.26675982 --output output\plot_adur2.json
python cli.py generate --lat 50.82989239 --lon -0.26747405 --output output\plot_adur3.json
python cli.py generate --lat 50.82992744 --lon -0.26734724 --output output\plot_adur4.json
python cli.py generate --lat 50.83296852141843 --lon -0.2758218218802195 --output output\plot_adur5.json
python cli.py generate --lat 51.268523 --lon -0.570246 --output output\plot_gu47nn2.json
python cli.py generate --lat 51.2678554 --lon -0.5695425105 --output output\plot_gu47nn3.json

python cli.py generate --lat 53.039271608492534 --lon -1.2006180639806219 --output output\plot_ash1.json


- Visualization Tests
python cli.py visualize --input output/plot_ash1.json --output output/plot_ash1.png



- Mapbox Imagery Tests
	# Basic usage (50m radius default)
python tests/test_mapbox_imagery.py --lat 51.268535 --lon -0.570979
	# Custom radius
python tests/test_mapbox_imagery.py --lat 51.26854 --lon -0.57042 --radius 40
	# Verbose mode for debugging
python tests/test_mapbox_imagery.py --lat 51.268535 --lon -0.570979 --verbose


- SAM3 segmentation Test
	# Download imagery and run SAM3
python tests/test_sam3_segmentation.py --lat 51.268535 --lon -0.570979
	# Custom radius
python tests/test_sam3_segmentation.py --lat 51.26854 --lon -0.57042 --radius 50
	# Use existing image (skip download)
python tests/test_sam3_segmentation.py --lat 51.268535 --lon -0.570979 --image tests/output/mapbox_merged_*.jpg
	# Verbose mode
python tests/test_sam3_segmentation.py --lat 51.268535 --lon -0.570979 --verbose


- Boundary Collector Test
	# Basic test
python tests/test_boundary_collector.py --lat 50.82991 --lon -0.266998
	# With custom radius and verbose output
python tests/test_boundary_collector.py --lat 51.268535 --lon -0.570979 --radius 50 --verbose
	# Test without OSM data pre-fetching
python tests/test_boundary_collector.py --lat 51.268535 --lon -0.570979 --no-osm-data


- Inspire boundary finder
python scripts/inspire_boundary_finder.py --lat 50.82991 --lon -0.266998
python scripts/inspire_boundary_finder.py --lat 50.82991 --lon -0.266998 --gml data/inspires/Adur_District_Council.gml
python scripts/inspire_boundary_finder.py --lat 50.82991 --lon -0.266998 --buffer 10 --output boundary.json
