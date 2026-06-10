# VASP NEB Monitor

A Python script to generate an HTML file to monitor VASP NEB (Nudged Elastic Band) calculation results.

## Overview

This tool helps you monitor and visualize the progress of VASP NEB calculations, which are used to find minimum energy paths and activation barriers between two atomic configurations.

## Features

- Generate interactive HTML monitoring reports
- Parse VASP NEB calculation results
- Clear and easy-to-read output format
- Cross-platform support (Windows, Linux)

## Usage

```bash
python neb_html_monitor.py <NEB_working_directory> [-o output_directory]
```

### Arguments

- `<NEB_working_directory>`: Path to your NEB working directory
- `-o, --output`: (Optional) Output file path for the HTML report

### Examples

**Windows:**
```bash
python neb_html_monitor.py F:\mace\prerun_gasS2_O2_1200K\30molecules\S_in
```

**Linux (current directory):**
```bash
python neb_html_monitor.py ./
```

**Linux (with custom output):**
```bash
python neb_html_monitor.py ./ -o ~/results.html
```

## Expected Directory Structure

Your NEB working directory should have the following structure:

```
NEB_working_directory/
├── INCAR
├── KPOINTS
├── POTCAR
├── 00/
│   ├── POSCAR
│   ├── CONTCAR
│   ├── OUTCAR
│   └── OSZICAR
├── 01/
│   ├── POSCAR
│   ├── CONTCAR
│   ├── OUTCAR
│   └── OSZICAR
└── ... (more image directories)
```

## Output

The script generates an HTML file that can be opened in any web browser. The output provides a clear and easy-to-read visualization of your VASP NEB calculation results.

## Requirements

- Python 3.x.x

## Installation

Clone the repository:

```bash
git clone https://github.com/wudiazhu/vasp-neb-monitor.git
cd vasp-neb-monitor
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

For questions or issues, please open an issue on GitHub.
