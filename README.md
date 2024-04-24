# UNAGI NEMO

GitHub repository of an idealised model I made in NEMO. See [here](https://nemo-related.readthedocs.io/en/latest/nemo_notes/unagi_config.html) for a description. The actual domain and forcing files may be found at the [zenodo repository](http://dx.doi.org/10.5281/zenodo.8002828). To save on uploading data however, it is possible to regenerate the files yourself through the provided notebooks; see also [here](https://nemo-related.readthedocs.io/en/latest/nemo_notes/unagi_config.html) for some instructions on how that might work (and is likely more transferable for making other NEMO models).

The model is used for the [Mak et al (2023) paper](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2023MS003886). See also the [zenodo repository](http://dx.doi.org/10.5281/zenodo.8002828) for a similar set of code but with the _splitting_ algorithm included.

# key notes and quick instructions:
* The GEOMETRIC code (see [personal repository](https://github.com/julianmak/GEOMETRIC_code) or the [NEMO Gitlab](https://forge.nemo-ocean.eu/nemo/nemo/-/tree/main/src/OCE/LDF?ref_type=heads)) mostly in `ldfeke.f90`.
* Wrote in by hand an enhanced vertical diffusivity profile for use with the channel model in `zdfphy.f90`.
* (BUG?) If model has linear ssh (which UNAGI does), it might default into the non-default computation for mean-flow advection of energy (which I personally don't remember and no longer understand...) The signature is that `trd_eke_adv_ubt` is smaller than all the other trends by about three orders of magnitude, and `eke` is too concentrated over the ridge. Code provided has the relevant parts commented out as a hack.
* (problem with TRUNK?) nemo4.2 trunk UNAGI for R100 seems to have a shift in the sponge region where vertical diffusivity is enhanced. Should have at least three points, which is the case in nemo4.2.1 and nemo4.2.2, but only one point in nemo4.2 trunk for some reason, using the same domcfg files.

# key updates:
* 24 Apr 2024 -- re-organised folders for (e.g. vanilla, with GEOMETRIC, with splitting) to work with different NEMO versions; includes the R100 configuration INPUT files (specifically for nemo4.2 domcfg files, which seem to work slightly differently)
* 22 Apr 2024 -- started working on the nemo4.2 versions (which seems to want a remake of the domcfg files; the forcing/state files seems ok)
* 08 Jun 2023 -- created repository
