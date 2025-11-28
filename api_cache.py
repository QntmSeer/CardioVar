import sqlite3
import json
import os
from datetime import datetime, timedelta
from typing import Optional, Any
import logging

class APICache:
    """Smart caching layer for API responses with TTL (Time To Live)."""
    
    def __init__(self, db_path: str = "data/api_cache.db", default_ttl_hours: int = 48):
        """
        Initialize the API cache.
        
        Args:
            db_path: Path to SQLite database file
            default_ttl_hours: Default cache lifetime in hours
        """
        self.db_path = db_path
        self.default_ttl_hours = default_ttl_hours
        
        # Ensure data directory exists
        os.makedirs(os.path.dirname(db_path), exist_ok=True)
        
        # Create table if it doesn't exist
        self._create_table()
    
    def _create_table(self):
        """Create the cache table if it doesn't exist."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS api_cache (
                cache_key TEXT PRIMARY KEY,
                data TEXT NOT NULL,
                created_at TIMESTAMP NOT NULL,
                expires_at TIMESTAMP NOT NULL
            )
        """)
        conn.commit()
        conn.close()
    
    def get(self, key: str, max_age_hours: Optional[int] = None) -> Optional[Any]:
        """
        Retrieve cached data if it exists and is not expired.
        
        Args:
            key: Cache key
            max_age_hours: Maximum age in hours (uses default if None)
        
        Returns:
            Cached data if valid, None otherwise
        """
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute("""
            SELECT data, expires_at FROM api_cache 
            WHERE cache_key = ?
        """, (key,))
        
        result = cursor.fetchone()
        conn.close()
        
        if not result:
            return None
        
        data_json, expires_at_str = result
        expires_at = datetime.fromisoformat(expires_at_str)
        
        # Check if expired
        if datetime.now() > expires_at:
            logging.debug(f"Cache expired for key: {key}")
            return None
        
        logging.debug(f"Cache hit for key: {key}")
        return json.loads(data_json)
    
    def set(self, key: str, data: Any, ttl_hours: Optional[int] = None):
        """
        Store data in cache with TTL.
        
        Args:
            key: Cache key
            data: Data to cache (must be JSON-serializable)
            ttl_hours: Time to live in hours (uses default if None)
        """
        ttl = ttl_hours if ttl_hours is not None else self.default_ttl_hours
        
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        created_at = datetime.now()
        expires_at = created_at + timedelta(hours=ttl)
        
        cursor.execute("""
            INSERT OR REPLACE INTO api_cache (cache_key, data, created_at, expires_at)
            VALUES (?, ?, ?, ?)
        """, (key, json.dumps(data), created_at.isoformat(), expires_at.isoformat()))
        
        conn.commit()
        conn.close()
        
        logging.debug(f"Cached data for key: {key} (TTL: {ttl}h)")
    
    def clear_expired(self):
        """Remove all expired entries from the cache."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute("""
            DELETE FROM api_cache 
            WHERE expires_at < ?
        """, (datetime.now().isoformat(),))
        
        deleted = cursor.rowcount
        conn.commit()
        conn.close()
        
        logging.info(f"Cleared {deleted} expired cache entries")
        return deleted
    
    def clear_all(self):
        """Clear all cache entries."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute("DELETE FROM api_cache")
        conn.commit()
        conn.close()
        logging.info("Cleared all cache entries")

    def invalidate(self, key: str):
        """Invalidate a specific cache key."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute("DELETE FROM api_cache WHERE cache_key = ?", (key,))
        conn.commit()
        conn.close()
        logging.info(f"Invalidated cache key: {key}")

    def invalidate_pattern(self, pattern: str):
        """Invalidate cache keys matching a SQL LIKE pattern."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute("DELETE FROM api_cache WHERE cache_key LIKE ?", (pattern,))
        deleted = cursor.rowcount
        conn.commit()
        conn.close()
        logging.info(f"Invalidated {deleted} keys matching pattern: {pattern}")
