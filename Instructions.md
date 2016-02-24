0. In `mint/` do `mkdir mint/projects`.
1. Put the `${project}_annotation.txt` file in `mint/projects/`. This informs the pipeline of the types of files and their relation to each other.
2. In `mint/`, run `bash init.sh project_name genome path_to_data` (creates project directory structure, moves the annotation file to `mint/projects/${project}/data/`, determines file names for `make` project runs, and symlinks the data from `path_to_data` into `mint/projects/${project}/data/raw_fastqs/`).
3. `cd projects/${project}`
4. In `mint/projects/${project}`, modify `config.mk` to reflect the project name and genome. Also modify the location of genomic indices, the paths to tools, and command line parameters for the pipeline.
5. In `mint/projects/${project}`, run `make copy_fastqs` (puts symlinks for data in `mint/${project}/data/raw_fastqs/` in `mint/${project}/bisulfite/raw_fastqs` and `mint/${project}/pulldown/raw_fastqs`, as necessary, with human- and pipeline-readable filenames).
6. Run `make bisulfite_align` and `make pulldown_align`, etc. Use the `-j` flag to use more processors.
