0. In `mint/` do `mkdir mint/projects`.
1. Put the `${project}_annotation.txt` file in `mint/projects/`. This informs the pipeline of the types of files and their relation to each other.
2. Ensure your data for the project are in one folder containing the `.fastq.gz`s.
3. In `mint/`, run `Rscript init.R --project project_name --genome genome --datapath path_to_data`
4. In `mint/projects/${project}`, modify `config.mk` to reflect the command line options for the tools. NOTE: The project name and genome are automatically populated from `init.R` options.
5. Modify `bisulfite_align.q`, `pulldown_align.q`, etc. for job parameters and `qsub` them.
