1. Put the `${project}_annotation.txt` file in `mint/projects/`. This informs the pipeline of the types of files and their relation to each other.
2. In `mint/`, run `bash init.sh project_name genome` (creates project directory structure, moves the annotation file to `mint/projects/${project}/data/`, and determines file names for `make` project runs).
3. Copy the raw `.fastq.gz` files into `mint/projects/${project}/data/raw_fastqs/` (symlink or full copy)
4. `cd projects/${project}`
5. In `mint/projects/${project}`, modify `config.mk` to reflect the project name and genome. Also modify the location of genomic indices, and command line parameters for the tools used in the pipeline.
6. In `mint/projects/${project}`, run `make copy_fastqs` (puts symlinks for data in `mint/${project}/data/raw_fastqs/` in `mint/${project}/bisulfite/raw_fastqs` and `mint/${project}/pulldown/raw_fastqs`, as necessary, with human- and pipeline-readable filenames).
7. Run `make bisulfite_align` and `make pulldown_align`, etc. Use the `-j` flag to use more processors.
