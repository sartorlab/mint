include project.mk

# Initialize a project
# $(PROJECT) comes from config.mk
.PHONY : init
init :
	@echo Making directories
	mkdir $(PROJECT)
	mkdir $(PROJECT)/data/
	mkdir $(PROJECT)/data/raw_fastqs
	mkdir $(PROJECT)/{bisulfite,pulldown}
	mkdir $(PROJECT)/bisulfite/{raw_fastqs,raw_fastqcs,trim_fastqs,trim_fastqcs,bismark,methylsig_calls}
	mkdir $(PROJECT)/pulldown/{raw_fastqs,raw_fastqcs,trim_fastqs,trim_fastqcs,bowtie2_bams,pulldown_coverages,macs_peaks,pepr_peaks}
	mkdir $(PROJECT)/{classification_simple,classification_sample,classification_comparison}
	mkdir $(PROJECT)/summary
	mkdir $(PROJECT)/summary/{figures,tables,reports}
	mkdir $(PROJECT)/$(PROJECT)_hub
	mkdir $(PROJECT)/$(PROJECT)_hub/$(GENOME)
	@echo Copying makefile and config.mk, and making empty variables.mk
	cp template_makefile $(PROJECT)/makefile
	cp template_config.mk $(PROJECT)/config.mk
	@echo Copying data annotation file to project data folder
	cp $(PROJECT)_annotation.txt $(PROJECT)/data/
	@echo Populating file lists
	awk -f scripts/parse_data.awk $(PROJECT)/data/$(PROJECT)_annotation.txt > $(PROJECT)/variables.mk
