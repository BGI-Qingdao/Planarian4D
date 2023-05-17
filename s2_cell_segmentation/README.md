### cell segmentation

- cell boundary was infered from ssDNA images using [CellProfiler](https://cellprofiler.org/) with custom parameters.

```shell
python3 cellsegment.py -i /path/to/ssDNA_image.tif -o ./outdir -p out_prefix --cpi default.cppipe
```

- for complex tasks of segmentation due to crowded cell environments or out-of-focus regions, Laplacian algorithm was used

```shell
python3 blur_detect.py /path/to/ssDNA_image.tif out_prefix
```

- finally, you can use the following to generation cell&bin mixtured mask, but you should generate a tissue segmentation mask file using [imageJ](https://imagej.net/ij/index.html) before running.

```shell
python3 complete_cell.py \
    --image /path/to/ssDNA_image.tif \
    --cell /path/to/cellprofiler_mask.txt \
    --blur /path/to/blurValue.txt \
    --tissue /path/to/ssDNA_tissue_mask.txt
```

