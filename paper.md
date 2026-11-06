---
title: 'gammapy-plugin: Including IACT data in threeML'
tags:
  - python
  - astronomy
  - astrophysics
  - high-energy
  - iact
authors:
  - name: Tobias Preis
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: 1
affiliations: 
 - name: University of Innsbruck
   index: 1

date: 
bibliography: paper.bib
---

# Summary
gammapy-plugin provides a plugin for the 3ML framework which enables joint Multi-Mission
analysis of Very High Energy data via the state-of-the-art software in this field, gammapy.

# Statement of need
Often MWL analysis of astrophysical sources rely on fitting flux points derived for
completely different models.
While there are current efforts in the gammapy community to allow for MWL analysis their
efforts are based on providing the data of all the different instrument in one common
data format.
3ML tries to solve this problem differently, while it also has support common data formats
like OGIP.
The approach is a plugin based framework, where each instrument has its own plugin, which 
acts as an interface to the official instrument software.
The plugin mainly get's handed a model, which its tasks is then to evaluate it instrument
specific and return a loglikelihood value which 3ML then uses for fitting and sampling.
This approach allows to fully exploit an instrument's capabilities without being restricted
by constraints from a common data format.
The latter is getting more relevant with e.g. polarization data. 

Including data from IACT instruments like e.g. H.E.S.S., MAGIC, Veritas or the upcoming
CTAO in 


# State of the field



# Software design
The general software design follows the 3ML approach of a plugin:
A model describing the region of interest gets passed to the plugin which then takes
care of folding the model through the relevant instrument responses and calculates the 
loglikelihood value.

In the case of the gammapy-plugin those steps are fully performed by gammapy itself.
The main task of this plugin is to convert the astromodel model into a gammapy model
and updating the parameters in the latter when they get changed during the fitting process.

Joint multiwavelength analysis in astrophysics plays a crucial role in constraining and
determining e.g. different emission mechanisms.
As different models can often reproduce similar observations, systematics from e.g. using
flux points produced for a different model, as little uncertainties as possible are 
required to determine the best fitting model.

While there exist solutions using a Common Data Format like OGIP to included other observations
directly in gammapy, 3ML uses the approach of providing a common interface to the official instrument
softwares.
This enables the user to not make any sacrifices 


# Research impact statement
This will plugin allows e.g. for joint analysis of IACT data with the native HAL plugin 
as well as the unbinned likelihood analysis for Fermi/LAT natively.

# Mathematics

# AI usage disclosure

Generative AI (namely GitHub Copilot and) was used for generating documentation,
optimizing existing code and generating new features. No generative AI was used for 
writing this manuscript.
For the code the authors of the commits containing AI generated code take full 
responsibility on functionality and security.

# Acknowledgements

# References
