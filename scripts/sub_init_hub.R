################################################################################
# Track hub files
# project_hub/genomes.txt
# project_hub/hub.txt
# project_hub/project_name_hub.html

# The trackDb.txt file is used no matter the situation
# Make sure we start fresh each time
hubtrackdbfile = sprintf('%s/trackDb.txt', dir_track, genome)
cat('', file=hubtrackdbfile)

# genomes.txt
genomes = c(
	sprintf('genome %s', genome),
	sprintf('trackDb %s/trackDb.txt', genome))
cat(genomes, file=sprintf('%s/genomes.txt', dir_hub), sep='\n')
# hub.txt
hub = c(
	sprintf('hub %s_hub', project),
	sprintf('shortLabel %s', project),
	sprintf('longLabel %s', project),
	'genomesFile genomes.txt',
	'email rcavalca@umich.edu',
	sprintf('descriptionUrl %s_hub.html', project))
cat(hub, file=sprintf('%s/hub.txt', dir_hub), sep='\n')
# project_name_hub.html
# TODO: Document what's in the hub!
	cat(' ', file=sprintf('%s/%s_hub.html',dir_hub, project))

priority = 1
# superTracks
if(bool_bis_comp || bool_pull_comp) {
	for(comparison in unique(comparisons$comparison)) {
		trackEntry = c(
			sprintf('track %s_group_comparison', comparison),
			'superTrack on show',
			sprintf('group %s_group_comparison', comparison),
			sprintf('shortLabel %s comparison', comparison),
			sprintf('longLabel %s comparison', comparison),
			sprintf('priority %s', priority),
			' ')
		cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
		priority = priority + 1
	}
}

if(bool_bis_samp || bool_pull_samp) {
	for(sample in sort(unique(samples$humanID))) {
		trackEntry = c(
			sprintf('track %s_sample', sample),
			'superTrack on show',
			sprintf('group %s_sample', sample),
			sprintf('shortLabel %s_sample', sample),
			sprintf('longLabel Sample-level tracks for %s', sample),
			sprintf('priority %s', priority),
			' ')
		cat(trackEntry, file=hubtrackdbfile, sep='\n', append=T)
		priority = priority + 1
	}
}
