# Plotting script for quality metrics of Emu abundance estimatations based on TRANA pipeline output

## Dependencies

- [Pixi](https://pixi.sh)
- [TRANA](https://github.com/genomic-medicine-sweden/TRANA)

## Installation

```bash
git clone https://github.com/kclinmicro/trana-abundance-qc-plotting.git
cd trana-abundance-qc-plotting
pixi init
```

## Usage

```
pixi shell
python main.py --input-dir <path-to-trana-output-dir>
```

Results will be created as .png plots alongside the `.tsv` output files from
TRANA/Emu.
