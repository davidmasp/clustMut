# CMC FORMAT

* position ( 1 based ) - same as GR and UCSC browser
* sample
* chromosome (problems with UCSC bla bla)
* alt
* ref (I think this could work as a validation and for "strand")

should somehow be divided by a separator.
This is a problem for speed.

```txt
chr1:48576549-48576549_T-C_TCGA-sample-barcode
```

It has 3 info fields,

position (in almost the bed format). The thing is that using the 1 based approach is super easy to load into the browser. Also from here you can get the MS96 or pentanucleotides..
substitution (it inherits the strand and the alternative base)
sample (then easy to mix different files)

split("_")

split("-")