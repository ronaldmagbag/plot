import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from src.collectors import SoilCollector

# Create soil collector instance
collector = SoilCollector()

# Test with UK coordinates
result = collector.get_soil_data(
    lat=53.025325345325484,
    lon=-1.2020914785758279
)

print(result)
