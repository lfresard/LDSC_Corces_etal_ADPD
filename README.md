# LD score regression: Cell type specific analyses

## 1. Change reference used for peaks to hg19
you can obtain the hg38ToHg19.over.chain.gz assembly conversion data needed [here](http://hgdownload.cse.ucsc.edu/goldenpath/hg38/liftOver/hg38ToHg19.over.chain.gz).
```
bash liftover_singlecell.sh <peak_directory>
```
