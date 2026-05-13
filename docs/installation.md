# Installation

You will need to install `threeML`, `astromodels` and of course `gammapy` within one
Python environment.
Thanks to `conda` this is now pretty straight forward.

## Conda/Mamba
```bash
conda install -c threeML -c conda-forge gammapy-plugin
```

Important here is the channel priority as the `conda-forge` `threeML` and `astromodels`
packages are not officially suppoorted by the threeML Team. Make sure those packages
come from the `threeML` channel.


## PyPi
`gammapy-plugin` is also available on PyPi and installable via pip. This should
take care of all Python dependencies.

```bash
pip install gammapy-plugin
```

## Source
If you read this you probably know how this goes:

```bash
git clone https://github.com/threeML/gammapy-plugin.git
cd gammapy-plugin
pip install .
```
