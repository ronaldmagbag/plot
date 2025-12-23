"""
Overpass API client

Handles communication with Overpass API including:
- Rate limiting
- Retry logic
- Error handling
"""

import time
import requests
from typing import Dict, Any
from loguru import logger

from ...config import get_config


class OverpassAPIClient:
    """Client for interacting with Overpass API"""
    
    def __init__(self):
        self.config = get_config()
        self.overpass_url = self.config.api.overpass_url
        self.timeout = self.config.api.overpass_timeout
        self._last_request_time = 0
        self._min_request_interval = 2.0  # Increased delay between requests
    
    def _rate_limit(self):
        """Ensure we don't exceed rate limits"""
        elapsed = time.time() - self._last_request_time
        if elapsed < self._min_request_interval:
            time.sleep(self._min_request_interval - elapsed)
        self._last_request_time = time.time()
    
    def query(self, query: str, retry_delay: float = 5.0) -> Dict[str, Any]:
        """
        Execute Overpass API query with retry logic
        
        Args:
            query: Overpass QL query string
            retry_delay: Initial delay between retries (increases with attempts)
            
        Returns:
            JSON response from Overpass API
            
        Raises:
            RuntimeError: If query fails after all retries
        """
        self._rate_limit()
        
        headers = {
            "User-Agent": self.config.api.user_agent,
            "Content-Type": "application/x-www-form-urlencoded"
        }
        
        for attempt in range(self.config.api.max_retries):
            try:
                response = requests.post(
                    self.overpass_url,
                    data={"data": query},
                    headers=headers,
                    timeout=self.timeout
                )
                response.raise_for_status()
                return response.json()
            except requests.exceptions.Timeout:
                wait_time = retry_delay * (attempt + 1)
                logger.warning(f"Overpass timeout (attempt {attempt + 1}/{self.config.api.max_retries}). Retrying in {wait_time}s...")
                if attempt < self.config.api.max_retries - 1:
                    time.sleep(wait_time)
                else:
                    logger.error(f"OSM API failed: Overpass timeout after {self.config.api.max_retries} attempts")
                    raise RuntimeError(f"Overpass API timeout after {self.config.api.max_retries} attempts")
            except requests.exceptions.HTTPError as e:
                if e.response.status_code in [429, 504] and attempt < self.config.api.max_retries - 1:
                    wait_time = retry_delay * (attempt + 1)
                    logger.warning(f"Overpass {e.response.status_code} (attempt {attempt + 1}/{self.config.api.max_retries}). Retrying in {wait_time}s...")
                    time.sleep(wait_time)
                else:
                    logger.error(f"OSM API failed: HTTP {e.response.status_code} after {self.config.api.max_retries} attempts")
                    raise RuntimeError(f"Overpass API HTTP error {e.response.status_code} after {self.config.api.max_retries} attempts") from e
            except requests.exceptions.RequestException as e:
                logger.warning(f"Overpass request failed (attempt {attempt + 1}): {e}")
                if attempt < self.config.api.max_retries - 1:
                    time.sleep(retry_delay * (attempt + 1))
                else:
                    logger.error(f"OSM API failed: Request exception after {self.config.api.max_retries} attempts: {e}")
                    raise RuntimeError(f"Overpass API request failed after {self.config.api.max_retries} attempts: {e}") from e
        
        return {"elements": []}

