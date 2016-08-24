function osl_recompile_fieldtrip()
%
% There seem to be recurrent problems caused by Mex files in FieldTrip.
% This script cleans up, compiles and replaces all occurrences of Mex files in FieldTrip.
%
% For additional issues involving mxSerialise, see: 
% http://undocumentedmatlab.com/blog/serializing-deserializing-matlab-data
%
% TL;DR: 
% Rename the extension to .cpp and replace prototypes declaration with:
%
% // MX_API_VER has unfortunately not changed between R2013b and R2014a,
% // so we use the new MATRIX_DLL_EXPORT_SYM as an ugly hack instead
% #if defined(__cplusplus) && defined(MATRIX_DLL_EXPORT_SYM)
%     namespace matrix{ namespace detail{ namespace noninlined{ namespace mx_array_api{
% #endif
% 
% EXTERN_C mxArray* mxSerialize(mxArray const *);
% EXTERN_C mxArray* mxDeserialize(const void *, size_t);
% // and so on, for any other MEX C functions that migrated to C++ in R2014a
% 
% #if defined(__cplusplus) && defined(MATRIX_DLL_EXPORT_SYM)
%     }}}}
%     using namespace matrix::detail::noninlined::mx_array_api;
% #endif
%
% JH

    global OSLDIR

    current_dir = pwd;
    fieldtrip_dir = fullfile( OSLDIR, 'spm12/external/fieldtrip' );
    
    % delete all source Mex-files
    src_clean    = find_and_delete_mex(fullfile( fieldtrip_dir, 'src' ));
    fileio_clean = find_and_delete_mex(fullfile( fieldtrip_dir, 'fileio/@uint64' ));
    
    % extract function names to check for conflicts
    src_names    = cellfun( @basename, src_clean, 'UniformOutput', false );
    fileio_names = cellfun( @basename, fileio_clean, 'UniformOutput', false );
    assert( isempty(intersect(src_names,fileio_names)), 'Conflicts between Mex function names.' );
    
    % map all other Mex files
    mex_files = find_mex_files_recursive(fieldtrip_dir);
    mex_clean = clean_mex_files(mex_files);
    
    % fix the mxErrMsgTxt problem
    getopt_file = fullfile( fieldtrip_dir, 'src/ft_getopt.c' );
    dk.fs.puts( getopt_file, strrep(dk.fs.gets(getopt_file),'mxErrMsgTxt','mexErrMsgTxt'), true );
    
    % recompile all source files
    cd(fieldtrip_dir);
    ft_compile_mex(true);
    cd(current_dir);
    
    % copy newly compiled version instead
    for i = 1:numel(mex_clean) 
        
        file = mex_clean{i};
        name = basename(file);
        
        switch name
            case src_names
                delete([file '.mex*']);
                copyfile( fullfile(fieldtrip_dir,'src',[name '.' mexext]), fileparts(file) );
            case fileio_names
                delete([file '.mex*']);
                copyfile( fullfile(fieldtrip_dir,'fileio/@uint64',[name '.' mexext]), fileparts(file) );
            otherwise
                warning('Could not find a replacement for Mex-file: "%s"',mex_clean{i});
        end
    end
    
end

function b = basename(file)
    [~,b]=fileparts(file); 
end

function mex_files = find_mex_files(folder)
    
    mex_files = dir(fullfile( folder, '*.mex*' ));
    mex_files = cellfun( @(x) fullfile(folder,x), {mex_files.name}, 'UniformOutput', false );
    
end

function mex_files = find_mex_files_recursive(folder)

    [~,mex_files] = system(['find "' folder '" -type f -name "*.mex*"']);
    mex_files = strsplit(mex_files,'\n');
    mex_files = mex_files(~cellfun( @isempty, mex_files, 'UniformOutput', true ));

end

function mex_map = find_and_delete_mex(folder)

    mex_files = find_mex_files(folder);
    mex_map   = clean_mex_files(mex_files);
    
    if ~isempty(mex_files)
        dk.println( '%d Mex-file(s) will be deleted:\n%s', numel(mex_files), strjoin(mex_files,'\n') );
        for i = 1:numel(mex_files)
            delete(mex_files{i});
        end
    else
        dk.println( 'No Mex-file found in folder "%s".', folder );
    end

end

function mex_clean = clean_mex_files(mex_files)

    n = numel(mex_files);
    mex_clean = cell(1,n);
    
    for i = 1:n
        [p,n,e] = fileparts(mex_files{i});
        mex_clean{i} = fullfile(p,n);
    end
    mex_clean = unique(mex_clean);
    
end
