Welcome to Gammapy Plugin's documentation!
===========================================================
.. image:: _static/banner-light.svg
   :class: only-light banner
   :alt: Project banner (light)

.. image:: _static/banner-dark.svg
   :class: only-dark banner
   :alt: Project banner (dark)


.. warning::
   When using `gammapy` MapDatasets and an `astromodels` ExtendedSource the resulting
   normalization of the spectrum will be a factor :math:`(180/\pi)^2` too small.
   This will be fixed within upcoming versions of `astromodels` and `threeML`.

.. warning::
   Currently `astromodels` only fully supports ICRS geometries.
   Please transform all your `gammapy` datasets accordingly.
   We are working on implementing this in `astromodels`.
   

About gammapy_plugin
--------------------


`Gammapy`_ Plugin is a plugin for `threeML`_ that allows to incorporate all
gammapy-supported instruments in your threeML analysis.

.. _Gammapy: https://gammapy.org/
.. _threeML: https://threeml.readthedocs.io/en/stable/


.. note::
   This project is still under very active development!

To get started check out the `get-started`_ section.

.. _get-started: intro.rst

.. toctree::
   :maxdepth: 1
   :caption: Contents
        
   intro
   notebooks/converter.ipynb

   api/API.rst

.. nblinkgallery::

   notebooks/crab_spectrum.ipynb
   notebooks/converter.ipynb
     
