# UNAGI NEMO

GitHub repository of an idealised model I made in NEMO. See [here](https://nemo-related.readthedocs.io/en/latest/nemo_notes/unagi_config.html) for a description. The actual domain and forcing files may be found at the [zenodo repository](http://dx.doi.org/10.5281/zenodo.8002828). To save on uploading data however, it is possible to regenerate the files yourself through the provided notebooks; see also [here](https://nemo-related.readthedocs.io/en/latest/nemo_notes/unagi_config.html) for some instructions on how that might work (and is likely more transferable for making other NEMO models).

The model is used for the paper TO ADD. See also the [zenodo repository](http://dx.doi.org/10.5281/zenodo.8002828) for a similar set of code but with the _splitting_ algorithm included.

# key notes and quick instructions:
* the GEOMETRIC code (see [personal repository](https://github.com/julianmak/GEOMETRIC_code) or the [NEMO Gitlab](https://forge.nemo-ocean.eu/nemo/nemo/-/merge_requests/202)) mostly in `ldfeke.f90`
* wrote in by hand an enhanced vertical diffusivity profile for use with the channel model in `zdfphy.f90`

# key updates:
* 08 Jun 2023 -- created repository
