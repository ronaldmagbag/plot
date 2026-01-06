
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


python cli.py generate --lat 53.040554321706324 --lon -1.2056779010008587 --output output\1.json
python cli.py generate --lat 53.04056394590633 --lon -1.2160624715683674 --output output\2.json
python cli.py generate --lat 53.04977822746581 --lon -1.2028364279935457 --output output\3.json

python cli.py generate --lat 53.0403338536591 --lon -1.2192989437839514 --output output\v1.json
python cli.py generate --lat 53.04085482610178 --lon -1.2190292066583432 --output output\v2.json
python cli.py generate --lat 53.03993081184168 --lon -1.219092465790886 --output output\v3.json
python cli.py generate --lat 53.03896051626191 --lon -1.2221818110142328 --output output\v4.json
python cli.py generate --lat 53.05236708155057 --lon -1.2026009414885899 --output output\v5.json

python cli.py generate --lat 50.9982535962166 --lon -0.10046479985243968 --output output\v6.json
python cli.py generate --lat 51.013456530550684 --lon -0.08674840435499474 --output output\plot.json
python cli.py generate --lat 51.01311213839883 --lon -0.08656023924511343 --output output\plot.json
python cli.py generate --lat 51.015173584845435 --lon -0.08716211166355772 --output output\plot.json
python cli.py generate --lat 51.01502404603471 --lon -0.0864140785181765 --output output\plot.json




51.01319565709778, -0.10537015106334707
python cli.py generate --lat 51.01319565709778 --lon -0.10537015106334707 --output output\plot.json

- Visualization Tests
python cli.py visualize --input output/plot.json --output output/plot.png
python scripts/draw_osm_on_mapbox.py --plot-json output/plot.json --cache-dir output/plot --output output/osm.jpg
python scripts/draw_property_on_mapbox.py --plot-json output/plot.json --cache-dir output/plot --output output/property.jpg


python .\tests\debug_classifier.py --input .\output\1.json --lat 53.040554321706324 --lon -1.2056779010008587 

python .\tests\debug_classifier.py --input .\output\v1.json --lat 53.0403338536591 --lon -1.2192989437839514 


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
python scripts/draw_property_on_mapbox.py --plot-json output/v1.json --cache-dir output/v1 --output output/v1_property.jpg

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
python scripts/draw_osm_on_mapbox.py --plot-json output/v1.json --cache-dir output/v1 --output output/v1_osm.jpg

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




