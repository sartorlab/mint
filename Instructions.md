0. In `mint/` do `mkdir mint/projects`.
1. Put the `${project}_annotation.txt` file in `mint/projects/`. This informs the pipeline of the types of files and their relation to each other.
2. In `mint/`, run `bash init.sh project_name genome path_to_data` (creates project directory structure, moves the annotation file to `mint/projects/${project}/data/`, determines file names for `make` project runs, and symlinks the data from `path_to_data` into `mint/projects/${project}/data/raw_fastqs/`).
3. `cd projects/${project}`
4. In `mint/projects/${project}`, modify `config.mk` to reflect the project name and genome. Also modify the location of genomic indices, the paths to tools, and command line parameters for the pipeline.
5. In `mint/projects/${project}`, run `make copy` (puts symlinks for data in `mint/${project}/data/raw_fastqs/` in `mint/${project}/bisulfite/raw_fastqs` and `mint/${project}/pulldown/raw_fastqs`, as necessary, with human- and pipeline-readable filenames).
6. In `mint/projects/${project}`, run `make pbs` to create job scripts for the bisulfite_align and pulldown_align parts. The number of processors, amount of ram, and email address for alerts need to be manually changed in `mint/projects/${project}/makefile` in the `pbs` rule.
7. `qsub bisulfite_align.q` and `qsub pulldown_align.q`
