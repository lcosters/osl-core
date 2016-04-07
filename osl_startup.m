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

if ~isdeployed

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check fsl has been setup
    assert( ~isempty('FSLDIR'), [ ...
        'The environmental variable FSLDIR is not set. Please exit Matlab, ensure FSLDIR is set and then restart Matlab. ' ...
        'See the Prerequisites section at https://sites.google.com/site/ohbaosl/download' ...
    ]);

    % try a dummy call to an fsl tool to make sure FSL is properly installed:
    [status,~] = system('fslval');
    assert( status ~= 0, [ ...
        'FSL is not installed properly. Perhaps check that the $FSLDIR/bin directory is in your PATH before starting Matlab. ' ...
        'See the Prerequisites section at https://sites.google.com/site/ohbaosl/download' ...
    ]);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save current path and restore it
    oldpath = strsplit(path,':');
    restoredefaultpath;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % check path for other versions of:
    checklist = { 'fieldtrip', 'spm', 'osl', 'mne', 'netlab', 'fsl', 'fmt' };
    oldpath   = oldpath(cellfun( @isempty, strfind(oldpath,matlabroot) )); % remove Matlab paths from oldpaths
    
    nold   = numel(oldpath);
    ncheck = numel(checklist);
    found  = false(1,nold);
    
    for i = 1:ncheck
        found = found | ~cellfun( @isempty, strfind(oldpath,checklist{i}) );
    end
    addpath(oldpath{~found});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % add external libraries
    external = { fullfile(getenv('FSLDIR'),'etc/matlab'), 'fmt', 'hmmbox_4_1', 'layouts', 'netlab3.3/netlab',...
        'netlab3.3/nethelp', 'ICA_tools/FastICA_25', 'ICA_tools/icasso122', 'spm12' };
    externalgen = {'std_masks'};
    
    cellfun( @(x) addpath(fullfile(osldir,x)), external );
    cellfun( @(x) addpath(genpath(fullfile(osldir,x))), externalgen );
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OSLDIR=osldir;

% Add osl2 at the end to override any other path (make sure we dont add the git folders though)
osl2_subfolders = strsplit(genpath(fullfile(osldir,'osl2')),':');
osl2_gitfolders = strfind( osl2_subfolders, '.git' );

addpath(osl2_subfolders{cellfun( @isempty, osl2_gitfolders )});

% Startup SPM12
spm_get_defaults('cmdline',true);
spm eeg;
close all;

% Show warning about removed paths at the end
if ~isdeployed && any(found)
    wmsg = strjoin( oldpath(found), '\n' );
    warning(['The following paths were removed because of conflicts:' 10 wmsg]);
end

% Remove fieldtrip replication
fpduplicate = fullfile( osldir, 'spm12/external/fieldtrip/external/' );
if exist(fpduplicate,'dir'), rmpath(genpath(fpduplicate)); end
