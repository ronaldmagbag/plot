
53.03931221785198, -1.200618125451693
53.03936634643617, -1.2000999470466172
53.039579267660336, -1.200269987188616
53.039407414357704, -1.2006222326947877


53.02564887164458, -1.201323337676318
53.029386664827584, -1.212482268566838
53.0256839932681, -1.2023543953647096
53.02546447440407, -1.2025623911542014
53.025325345325484, -1.2020914785758279

- Generation Plot Tests
python cli.py generate --lat 50.82991 --lon -0.266998 --output output\plot_adur1.json
python cli.py generate --lat 50.82981237 --lon -0.26675982 --output output\plot_adur2.json
python cli.py generate --lat 50.82989239 --lon -0.26747405 --output output\plot_adur3.json
python cli.py generate --lat 50.82992744 --lon -0.26734724 --output output\plot_adur4.json
python cli.py generate --lat 50.83296852141843 --lon -0.2758218218802195 --output output\plot_adur5.json
python cli.py generate --lat 51.268523 --lon -0.570246 --output output\plot_gu47nn2.json
python cli.py generate --lat 51.2678554 --lon -0.5695425105 --output output\plot_gu47nn3.json



python cli.py generate --lat 53.03931221785198 --lon -1.200618125451693 --output output\plot_ash1.json
python cli.py generate --lat 53.03936634643617 --lon -1.2000999470466172 --output output\plot_ash2.json
python cli.py generate --lat 53.039579267660336 --lon -1.200269987188616 --output output\plot_ash3.json


python cli.py generate --lat 53.02564887164458 --lon -1.201323337676318 --output output\plot_ash4.json
python cli.py generate --lat 53.029386664827584 --lon -1.212482268566838 --output output\plot_ash5.json
python cli.py generate --lat 53.0256839932681 --lon -1.2023543953647096 --output output\53.0256839932681-1.2023543953647096.json
python cli.py generate --lat 53.025325345325484 --lon -1.2020914785758279 --output output\53.025325345325484-1.2020914785758279.json



- Visualization Tests
python cli.py visualize --input output/53.025325345325484-1.2020914785758279.json --output output/53.025325345325484-1.2020914785758279.png


python .\tests\debug_classifier.py --input .\output\53.025325345325484-1.2020914785758279.json --lat 53.025325345325484 --lon -1.2020914785758279 


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


- draw_property_on_mapbox
# Basic usage - auto-find image from cache
python scripts/draw_property_on_mapbox.py --plot-json output/plot_ash1.json --cache-dir output/plot_ash1 --output output/plot_ash1_with_property.jpg

# Specify exact image
python scripts/draw_property_on_mapbox.py \
    --image output/plot_ash1/mapbox_merged_53.039271_-1.200618_20.0m_z19.jpg \
    --plot-json output/plot_ash1.json \
    --output output/annotated.jpg

# Custom styling
python scripts/draw_property_on_mapbox.py \
    --plot-json output/plot_ash1.json \
    --cache-dir output/plot_ash1 \
    --output output/annotated.png \
    --line-color "0,255,0" \
    --line-width 5 \
    --image-alpha 0.6



- Draw all OSM features with default colors
python scripts/draw_osm_on_mapbox.py --plot-json output/plot_ash1.json --cache-dir output/plot_ash1 --output output/plot_ash1_osm.jpg

# Draw only buildings and roads
python scripts/draw_osm_on_mapbox.py \
    --plot-json output/plot_ash1.json \
    --cache-dir output/plot_ash1 \
    --output output/plot_ash1_osm.jpg \
    --no-trees --no-water

# Custom colors and styling
python scripts/draw_osm_on_mapbox.py \
    --plot-json output/plot_ash1.json \
    --cache-dir output/plot_ash1 \
    --output output/plot_ash1_osm.jpg \
    --building-color "100,100,100" \
    --road-color "255,255,0" \
    --tree-color "0,200,0" \
    --road-width 5 \
    --tree-radius 6 \
    --image-alpha 0.6

# Specify exact image file
python scripts/draw_osm_on_mapbox.py \
    --image output/plot_ash1/mapbox_merged_53.039271_-1.200618_20.0m_z19.jpg \
    --plot-json output/plot_ash1.json \
    --output output/annotated_osm.jpg