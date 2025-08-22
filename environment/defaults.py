from pathlib import Path

DEFAULT_CONTAINER_CONFIG_PATH = (
    Path(__file__).parent / "config" / "container_configs" / "default.json"
)

CPU_ONLY_CONTAINER_CONFIG_PATH = (
    Path(__file__).parent / "config" / "container_configs" / "cpu-only.json"
)

FAST_CONTAINER_CONFIG_PATH = (
    Path(__file__).parent / "config" / "container_configs" / "fast.json"
)
