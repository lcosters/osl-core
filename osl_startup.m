function osl_startup( osldir )

global OSLDIR;

% osl_startup( osldir )
%
% SETS UP THE BASIC PATH SETTINGS
%
% MWW 2012
%
% does no path-changing if running in deployed mode (gw '13).
% modernise the code (jh 2016)

    fprintf( '[OSL2] Starting up from folder "%s"...\n', osldir );

    if ~isdeployed

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf( '\t Check FSL installation...\n' );
        
        % check fsl has been setup
        assert( ~isempty('FSLDIR'), [ ...
            'The environment variable FSLDIR is not set. Please exit Matlab, ensure FSLDIR is set and then restart Matlab. ' ...
            'See the Prerequisites section at https://sites.google.com/site/ohbaosl/download' ...
        ]);

        % try a dummy call to an fsl tool to make sure FSL is properly installed:
        [status,~] = system('fslval');
        assert( status ~= 0, [ ...
            'FSL is not installed properly. Perhaps check that the $FSLDIR/bin directory is in your PATH before starting Matlab. ' ...
            'See the Prerequisites section at https://sites.google.com/site/ohbaosl/download' ...
        ]);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf( '\t Set paths...\n' );
        
        % set path
        inc_ext = { 'fmt', 'hmmbox_4_1', 'layouts', 'netlab3.3/netlab', 'netlab3.3/nethelp', ...
            'ICA_tools/FastICA_25', 'ICA_tools/icasso122', 'spm12' };
        inc_ext = cellfun( @(x) fullfile(osldir,x), inc_ext, 'UniformOutput', false );
        inc_ext = strjoin( inc_ext, pathsep );

        inc_osl = { 'africa', 'HCP', 'oat', 'oil', 'opt', 'osl_hmm_toolbox', 'oslview', 'rhino', 'util' };
        inc_osl = cellfun( @(x) genpath(fullfile( osldir, 'osl2', x )), inc_osl, 'UniformOutput', false );
        inc_osl = [ fullfile(osldir,'osl2'), pathsep, strjoin(inc_osl,pathsep) ];

        includes = { inc_osl, genpath(fullfile( osldir, 'std_masks' )), inc_ext };
        includes = strjoin( includes, pathsep );
        includes = strsplit( includes, pathsep );

        curpath = strsplit(path,pathsep);
        newpath = unique( [includes,curpath] );
        path(strjoin( newpath, pathsep ));

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf( '\t Apply changes to SPM...\n' );
    
    % Copy changes to SPM code from osl
    filelist={};targetdir={};

    filelist{end+1}='osl2/spm-beamforming-toolbox-osl-addons/bf_output_montage_osl.m';
    targetdir{end+1}='spm12/toolbox/spm-beamforming-toolbox';

    filelist{end+1}='osl2/spm-beamforming-toolbox-osl-addons/bf_write_spmeeg_osl.m';
    targetdir{end+1}='spm12/toolbox/spm-beamforming-toolbox';

    filelist{end+1}='osl2/spm-beamforming-toolbox-osl-addons/bf_inverse_mne_adaptive.m';
    targetdir{end+1}='spm12/toolbox/spm-beamforming-toolbox';

    filelist{end+1}='osl2/spm-changes/private/ft_read_event_4osl.m';
    targetdir{end+1}='spm12/external/fieldtrip/fileio';

    filelist{end+1}='osl2/spm-changes/private/read_trigger_4osl.m';
    targetdir{end+1}='spm12/external/fieldtrip/fileio/private';

    filelist{end+1}='osl2/spm-changes/private/badsamples.m';
    targetdir{end+1}='spm12/@meeg';

    filelist{end+1} = 'osl2/spm-changes/private/undobalancing.m';
    targetdir{end+1}='spm12/external/fieldtrip/forward/private/';

    filelist{end+1} = 'osl2/spm-changes/private/undobalancing.m';
    targetdir{end+1}='spm12/external/fieldtrip/plotting/private/';

    filelist{end+1} = 'osl2/spm-changes/private/undobalancing.m';
    targetdir{end+1}='spm12/external/fieldtrip/fileio/private/';

    filelist{end+1} = 'osl2/spm-changes/private/undobalancing.m';
    targetdir{end+1}='spm12/external/fieldtrip/utilities/private/';

    filelist{end+1} = 'osl2/spm-changes/private/undobalancing.m';
    targetdir{end+1}='spm12/external/fieldtrip/private/';

    filelist{end+1} = 'osl2/spm-changes/private/ft_headmodel_localspheres.m';
    targetdir{end+1}='spm12/external/fieldtrip/forward/';

    filelist{end+1}='osl2/spm-changes/private/path.m';
    targetdir{end+1}='spm12/@meeg';

    filelist{end+1}='osl2/spm-changes/private/spm_eeg_montage.m';
    targetdir{end+1}='spm12';

    filelist{end+1} ='osl2/spm-changes/private/spm_eeg_inv_mesh_ui.m';
    targetdir{end+1}='spm12';

    filelist{end+1} ='osl2/spm-changes/private/subsref.m';
    targetdir{end+1}='spm12/@meeg';

    filelist{end+1} ='osl2/spm-changes/private/ft_getopt.c';
    targetdir{end+1}='spm12/external/fieldtrip/src/';

    for k=1:length(filelist),
        copyfile( fullfile(osldir,filelist{k}), fullfile(osldir,targetdir{k}), 'f' );
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    OSLDIR=osldir;

    % Startup SPM12
    fprintf( '\t Start SPM...\n' );
    spm_jobman('initcfg');
    spm('ChMod','eeg');
    close all;

end
