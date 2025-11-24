from gammapy_plugin.utils.package_data import (
    get_path_of_data_dir,
    get_path_of_data_file,
)


def test_get_path_of_data_dir():
    data_dir = get_path_of_data_dir()
    assert data_dir.is_dir()
    assert get_path_of_data_dir() / "create_test_datasets.py" in list(
        get_path_of_data_dir().iterdir()
    )


def test_get_path_of_data_file():
    assert get_path_of_data_file("create_test_datasets.py").is_file()
