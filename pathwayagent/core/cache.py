"""Tiny SQLite cache for external-service responses (design §7.4)."""
from __future__ import annotations

import hashlib
import json
import sqlite3
import threading
import time

from ..config import CACHE_PATH

_lock = threading.Lock()
TTL_S = 14 * 24 * 3600


def _conn():
    CACHE_PATH.parent.mkdir(parents=True, exist_ok=True)
    c = sqlite3.connect(CACHE_PATH, timeout=30)
    c.execute("CREATE TABLE IF NOT EXISTS cache (k TEXT PRIMARY KEY, t REAL, v TEXT)")
    return c


def make_key(*parts) -> str:
    return hashlib.sha256(json.dumps(parts, sort_keys=True, default=str).encode()).hexdigest()


def get(key: str):
    with _lock:
        c = _conn()
        row = c.execute("SELECT t, v FROM cache WHERE k=?", (key,)).fetchone()
        c.close()
    if row and time.time() - row[0] < TTL_S:
        return json.loads(row[1])
    return None


def put(key: str, value) -> None:
    with _lock:
        c = _conn()
        c.execute("INSERT OR REPLACE INTO cache VALUES (?,?,?)", (key, time.time(), json.dumps(value)))
        c.commit()
        c.close()


def cached(namespace: str, payload, fn, enabled: bool = True, store_if=bool):
    """Return fn() result, cached by (namespace, payload). Empty/failed results are not stored."""
    key = make_key(namespace, payload)
    if enabled:
        hit = get(key)
        if hit is not None:
            return hit, True
    val = fn()
    if enabled and store_if(val):
        put(key, val)
    return val, False
