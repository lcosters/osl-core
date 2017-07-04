function osl_shutdown
	% Remove OSL from path
	%
	% Note that this will not remove FSL or workbench from the system path
	% so commands like "system('fslval')" should still run
	% In that respect, there will still be some traces of OSL. It's a bit more 
	% complicated to remove these because they may not have been added by
	% OSL, so best not to worry about them for now. The FSL Matlab helpers
	% may also be left intact, if they were present before OSL was started
	%
	% Romesh Abeysuriya 2017

	% OSL is initialized if OSLDIR is set. If OSLDIR is not set, then do nothing
	if isempty(getenv('OSLDIR'))
		return
	else
		osl_root = getenv('OSLDIR');
	end

	% If path backup exists, read it
	if exist(fullfile(osl_root,'path_backup.txt'),'file')
		fid = fopen(fullfile(osl_root,'path_backup.txt'),'r');
		backup_path = fgetl(fid);
		path(backup_path);
		fclose(fid);
		delete(fullfile(osl_root,'path_backup.txt')); % Delete the backup file now that the path has been reset
	end

	% It's concievable that the backup path still contains OSL directories, if osl_startup
	% was run twice. Remove them just in case
	p = regexp(path,':','split');
	for j = 1:length(p)
		if strfind(p{j},osl_root)
			rmpath(p{j});
		end
	end

	% Done - clear OSLDIR
	setenv('OSLDIR');
