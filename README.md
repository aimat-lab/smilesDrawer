# SMILES to image
Use [generate-images.sh](generate-images.sh) to generate data. 
You need to install [node.js](https://nodejs.org/en) to run this tool.

Install dependencies using
```console
npm i
```

| Parameter             | Description                                                                   |
|-----------------------|-------------------------------------------------------------------------------|
| `--from-csv-file`     | Path to a CSV file. It can have one column only.                              |
| `--from-csv-column`   | Specifies which column of CSV to read. Set to 0 for CVS with only one column. |
| `--output-directory`  | Output directory.                                                             |
| `--size`              | Value specifying by how much the generated image should be resized.           |
| `--fonts`             | Fonts to use. Can be any font supported by Chromium.                          |
| `--font-weights`      | Fonts weights to use.                                                         |
| `--concurrency`       | How many headless browsers to start.                                          |
| `--min-smiles-length` | Lower bound for SMILES strings.                                               |
| `--max-smiles-length` | Upper bound for SMILES strings.                                               |
| `--amount`            | How many SMILES to read.                                                      |
| `--batch-size`        | How many images to generate in parallel.                                      |
| `--output-labels`     | Whether to output labels.                                                     |
| `--output-svg`        | Debug option. Whether to output raw SVG files.                                |
| `--output-flat`       | Debug option. Whether to output all files into the same directory.            |
| `--clean`             | Debug option. Whether to clean the target directory.                          |



