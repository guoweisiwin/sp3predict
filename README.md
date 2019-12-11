# sp3predict

Simple script to make DST predictions given a (Clockwork) VCF, a catalogue and a gumpy genome object.


## Pre-requisities

1. `gumpy`.
2. `piezo`
3. One or more AMR catalogues in the format required by `piezo`. Several published tuberculosis catalogues are available (here)[https://github.com/oxfordmmm/tuberculosis_amr_catalogues]

Since the latter is large, rather than store it in the repo, recreate the object using `gumpy-save-genome.py` which should be in your `$PATH` if you have already installed `gumpy`.

```
$ gumpy-save-genome.py --genbank_file config/H37rV_v3.gbk --name config/H37rV_v3.pkl
```

This will take about five minutes.

Now you are ready 
