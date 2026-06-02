# A collection of tools and scripts for plotting and troubleshooting output from the Emu tool

## Dependencies

- [Pixi](https://pixi.sh)
- [Emu](https://github.com/treangenlab/emu)
- ([TRANA](https://github.com/genomic-medicine-sweden/TRANA))

## Installation

```bash
git clone https://github.com/kclinmicro/emutools.git
cd emutools
pixi init
```

## Usage: Plotting probability densities

```
pixi shell
python plot_probdist.py --input-dir <path-to-trana-output-dir>
```

Results will be created as .png plots alongside the `.tsv` output files from
TRANA/Emu.
