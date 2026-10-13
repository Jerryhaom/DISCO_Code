# Current submission analyses

These modules accompany Fig. 3c / Supplementary Table 16 and Supplementary Fig. S8 / Table 26E3. They retain the models verified on 11 October 2026. Earlier bundled modules remain as historical reproductions and do not supersede these updates.

* `dm`: natural-log DM, standardized in the full scored population of each modality; crude and eight-covariate disease HRs. The figure retains crude score-only C-indices.
* `matched`: DISCO on each clock-selected protein set, with seven-covariate and mutually adjusted mortality models; scores use frozen training age-residual SDs.

Participant records are not distributed. Authorized users need the corresponding Results_final input tree listed in the R scripts. `PROVENANCE.json` records the exact audited source hashes. Analytic expressions are unchanged; paths and plotting devices are portable.

```sh
python3 submission_20261013/run.py dm --data-root /path/to/Results_final --output /path/to/new-dm-run --check-inputs
python3 submission_20261013/run.py dm --data-root /path/to/Results_final --output /path/to/another-new-dm-run
python3 submission_20261013/run.py matched --data-root /path/to/Results_final --output /path/to/new-matched-run
```

Use R with data.table, survival, jsonlite and digest; plotting additionally requires forestploter, gtable, ggplot2 and patchwork. Optionally set DISCO_R_LIBS. Cairo support and Arial are needed to reproduce the figure typography. Input checks validate files and parse code; they do not refit models. Outputs are restricted to a new directory.
