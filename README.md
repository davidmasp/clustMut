# clustMut

The goal of clustMut is to ...

## Installation

Install the package using:

```bash
...
```

Then, move the script to run to your `bin` or to a folder available in your `$PATH`.

```bash
cp clustmut_run.sh ~/bin/
```

and give it execution permisions.

```bash
chmod +x script.sh
```

## Usage

You can use clustmut to obtain kataegis events, clusters based on VAF or clustered mutations based on the Edit distance.
Run it with the appropiate mode command.

### VAF

```bash
clustmut_run.sh --mode vaf \
                -i /home/dmas/data/TCGA_MUTS/TCGA_VR/ \
                --glob "*_VR.rds" \
                --recursive \
                -a ~/data/CRG_alignability/hg19/LEGACY/crg36AlignExtNoBlackRmvsh19_RngMask_savedInt\=TRUE.bed \
                -o test \
                -Vlwtvu
```

### distance (difuse clusters)

```bash
clustmut_run.sh --mode distance \
                -i /home/dmas/data/TCGA_MUTS/RNDmut/ \
                --glob "*randomized.tsv" \
                --recursive \
                -o test_omichili \
                -N 1 \
                -Vlwtvu
```

### distance (kataegis)

```bash
clustmut_run.sh --mode distance \
                -i /home/dmas/data/TCGA_MUTS/RNDmut/ \
                --glob "*randomized.tsv" \
                --recursive \
                -o test_kataegis \
                -N 4 \
                -Vlwtvu
```

### edit

```bash
clustmut_run.sh --mode edit \
                -i /home/dmas/data/TCGA_MUTS/TCGA_VR/ \
                --glob "*_VR.rds" \
                --recursive \
                -a ~/data/CRG_alignability/hg19/LEGACY/crg36AlignExtNoBlackRmvsh19_RngMask_savedInt\=TRUE.bed \
                -o test \
                -Vlwtvu
```
