# vasp-neb-monitor
A script to generate html file to monitor the VASP NEB results

Usage:
python neb_html_monitor.py <NEB_working_directory> -o [output_directory]

Example:
(WIN)   python neb_html_monitor.py F:\mace\prerun_gasS2_O2_1200K\30molecules\S_in
(LINUX) python neb_html_monitor.py ./
(LINUX) python neb_html_monitor.py ./ -o ~/results.html

Expected NEB Directory:
NEB_working_directory/
  INCAR
  KPOINTS
  POTCAR
  00/
    POSCAR
    CONTCAR
    OUTCAR
    OSZICAR
  01/
    POSCAR
    CONTCAR
    OUTCAR
    OSZICAR
  ...

See the output:
Open the html file with your browser. The output content is clearly defined and easy to read.
