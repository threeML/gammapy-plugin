from dataclasses import dataclass
from pathlib import Path

from omegaconf import OmegaConf
from omegaconf.dictconfig import DictConfig
from rich.tree import Tree

# Path to configuration

_config_path = Path("~/.config/gammapy_plugin/").expanduser()

_config_name = Path("gammapy_plugin_config.yml")

_config_file = _config_path / _config_name

# Define structure of configuration with dataclasses


@dataclass
class Logging:

    on: bool = True
    level: str = "WARNING"


@dataclass
class gammapy_pluginConfig:

    logging: Logging = Logging()


# Read the default config
gammapy_plugin_config: gammapy_pluginConfig = OmegaConf.structured(gammapy_pluginConfig)

# Merge with local config if it exists
if _config_file.is_file():

    _local_config = OmegaConf.load(_config_file)

    gammapy_plugin_config: gammapy_pluginConfig = OmegaConf.merge(
        gammapy_plugin_config, _local_config
    )

# Write defaults if not
else:

    # Make directory if needed
    _config_path.mkdir(parents=True, exist_ok=True)

    with _config_file.open("w") as f:

        OmegaConf.save(config=gammapy_plugin_config, f=f.name)


def recurse_dict(d, tree) -> None:

    for k, v in d.items():

        if (type(v) is dict) or isinstance(v, DictConfig):

            branch = tree.add(
                k, guide_style="bold medium_orchid", style="bold medium_orchid"
            )

            recurse_dict(v, branch)

        else:

            tree.add(
                f"{k}: [blink cornflower_blue]{v}",
                guide_style="medium_spring_green",
                style="medium_spring_green",
            )

    return


def show_configuration() -> None:

    tree = Tree("config", guide_style="bold medium_orchid", style="bold medium_orchid")

    recurse_dict(gammapy_plugin_config, tree)

    return tree
