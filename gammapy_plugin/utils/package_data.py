import importlib
from pathlib import Path

__all__ = [
    "get_path_of_data_dir",
    "get_path_of_data_file",
]


def get_path_of_data_dir() -> Path:
    """Get the path of the package data directory.

    :returns:
    """
    file_path = importlib.resources.files("gammapy_plugin.data")

    return file_path


def get_path_of_data_file(data_file: str) -> Path:
    """Get the path of a data file.

    :param data_file: name of the data file
    :type data_file: str
    :returns:
    """
    file_path: Path = get_path_of_data_dir() / data_file
    return file_path
