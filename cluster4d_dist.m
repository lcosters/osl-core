function stats = cluster4d_dist(S);

% subfunction of cluster4d_batch
% takes the output of cluster4d_batch and generates a null distribution
%   against which to compare cluster sizes from the true T-statistic, and
%   produces stats images
%
% input arguments:
% try, dirname = S.dirname; catch, error('Must specify S.dirname'); end
% try, tp = S.tp; catch, error('Must specify S.tp'); end
% try, np = S.np; catch, np = 1:5000; end
% try, distfile_save = S.distfile_save; catch, distfile_save = [];end
% try, distfile_load = S.distfile_load; catch, distfile_load = [];end
% try, save_images = S.save_images; catch, save_images =0; end
%
% Laurence Hunt and Tom Nichols and Mark Woolrich, 2010/11 (beta version!)


try, dirname = S.dirname; catch, error('Must specify S.dirname'); end
try, subdirname = S.subdirname; catch, error('Must specify S.subdirname'); end
try, tp = S.tp; catch, error('Must specify S.tp'); end
try, times = S.times; catch, error('Must specify S.times'); end
try, np = S.np; catch, np = 1:5000; end
try, distfile_save = S.distfile_save; catch, distfile_save = [];end
try, distfile_load = S.distfile_load; catch, distfile_load = [];end
try, save_images = S.save_images; catch, save_images =1; end
try, gridstep = S.gridstep; catch, gridstep=2;  end

nT =length(tp);

%% make clusters from true T-stat image
first = 1;
%load in the true tstat

tstats=[];
for i = 1:length(tp);
    if strcmp(S.permmeth, 'clustextent')
        fname = sprintf('%s/%s/%04.0f/stats%04.0f_clustere_tstat1.nii.gz',dirname,subdirname,tp(i),tp(i));
    elseif strcmp(S.permmeth, 'clustmass')
            fname = sprintf('%s/%s/%04.0f/stats%04.0f_clusterm_tstat1.nii.gz',dirname,subdirname,tp(i),tp(i));
    end
    f = nii.load(fname);
    fname_rawt = sprintf('%s/%s/%04.0f/stats%04.0f_tstat1.nii.gz',dirname,subdirname,tp(i),tp(i));
    tmp=nii.load(fname_rawt);
    if(isempty(tstats)),
       tstats=zeros([size(tmp) length(tp)]); 
    end;
    
    try 
        tstats(:,:,:,i) = tmp;
    catch 
        tstats(:,:,i)=tmp; %small workaround for working with parcellated data
    end

    if first %first image - read in dimensions
        Vdim = size(f);
        nV = prod(Vdim);
        Xreal = sparse(nV,nT); % XT is nVoxels x timepoints and contains t-values of clusters (needed for clustermass permutation tests)
        first = 0;
    end
    Xreal(:,i) = f(:);
end
Xreal = clusterX(Xreal); % Xreal is nVoxels x timepoints and contains cluster sizes (clustextent) or cluster t-values (clustmass)
Creal = unique(Xreal);
Creal(Creal==0) = [];
if ~isempty(Creal)
    for i = 1:length(Creal)
        nVreal(i) = sum(Xreal(find(Xreal))==Creal(i)); 
    end
else
    nVreal = 0;
end
clustimg = (reshape(full(Xreal),[Vdim nT])); % clustimg is 4D matrix containing the cluster sizes or t-values for all (merged) clusters

%% create null distribution
if isempty(distfile_load) %no distribution to load, so we create it from scratch
    if strcmp(S.permmeth, 'clustextent')
        dist = zeros(length(np),1);
        for p = 1:length(np) %loop over permutations
            
            first = 1;
            %% first, load in files
            for i = 1:length(tp);
                fname = sprintf('%s/%s/%04.0f/cindex%05.0f',dirname,subdirname,tp(i),np(p));
                f = nii.load(fname);
                if first %first image - read in dimensions
                    nV = prod(size(f));
                    Xtmp = sparse(nV,nT); % XT is nVoxels x timepoints and contains t-values of clusters (needed for clustermass permutation tests)
                    first = 0;
                end
                Xtmp(:,i) = f(:);
            end
            
            Xtmp = clusterX(Xtmp); % XI is nVoxels x timepoints and contains indexes of clusters (which are simply labels and have no meaning)
            Ctmp = unique(Xtmp);
            Ctmp(Ctmp==0) = [];
            if ~isempty(Ctmp)
                for i = 1:length(Ctmp)
                    nVtmp(i) = sum(Xtmp(find(Xtmp))==Ctmp(i)); % nV is nClusters x 1 and contains the number of voxels in each cluster (needed for cluster extent permutation tests)
                end
            else
                nVtmp = 0;
            end
            %now get the maximum number of voxels in a single cluster, and add to dist
            dist(p) = max(nVtmp);
            clear tmp;
            fprintf('Permutation %0.0f complete...\n',p);
        end
        
    elseif strcmp(S.permmeth, 'clustmass')
        dist = zeros(length(np),1);
        for p = 1:length(np) %loop over permutations
            first = 1;
            %% first, load in files
            for i = 1:length(tp);
                %fname = sprintf('%s/%s/%04.0f/cindex%05.0f',dirname,subdirname,tp(i),np(p));
                fname = sprintf('%s/%s/%04.0f/stats%04.0f_clusterm_tstat1_perm%05d.nii.gz',dirname,subdirname,tp(i),tp(i), np(p));
                f = nii.load(fname);
                if first %first image - read in dimensions
                    nV = prod(size(f));
                    Xtmp = sparse(nV,nT); % XT is nVoxels x timepoints and contains t-values of clusters (needed for clustermass permutation tests)
                    first = 0;
                end
                Xtmp(:,i) = f(:);
            end
            
            Xtmp = clusterX(Xtmp); % XI is nVoxels x timepoints and contains indexes of clusters (which are simply labels and have no meaning)
            Ctmp = unique(Xtmp);
            dist(p) = max(Ctmp);
            clear tmp;
            fprintf('Permutation %0.0f complete...\n',p);
        end
    end
else
    if ~iscell(distfile_load)
        error('S.distfile_load must be a cell array');
        %error('Not yet implemented loading of distributions');
    else
        for i = 1:length(distfile_load)
            dist = [];
            Ctmp = load(distfile_load{i});
            dist = [dist Ctmp.dist];
        end
    end
end

%% see where clusters lie in distribution, assign p-values and write images

if save_images %not yet implemented for cluster mass
    
    nC = length(nVreal);
    pVreal = zeros(nC,1);
    pVimg = clustimg;
    for i = 1:nC %loop over real clusters
        pVreal(i) = mean(nVreal(i)>=dist);
        pVimg(clustimg==full(Creal(i))) = full(pVreal(i));
    end

    clustimg_fname = [dirname '/clust4d_' num2str(min(np)) '-' num2str(max(np)) '.nii.gz']; 
    pVimg_fname = [dirname '/clust4d_corrp_' num2str(min(np)) '-' num2str(max(np)) '.nii.gz']; %make clear that the files that are saved are only for 100 permutations
    
    xform = [-gridstep 0 0 90; 0 gridstep 0 -126; 0 0 gridstep -72; 0 0 0 1]; 
    tres=1;
    nii.save(clustimg,[gridstep gridstep gridstep tres],xform,clustimg_fname);
    nii.save(pVimg,[gridstep gridstep gridstep tres],xform,pVimg_fname);
    
end

if ~isempty(distfile_save)
    save(distfile_save,'nVreal','clustimg','Creal','dist','times','tstats');
end
end

function Xin=clusterX(Xin);
tic
nV = size(Xin,1);
nT = size(Xin,2);
% Xin = nV by nT matrix - input matrix & directly modified to create output matrix
% On input, Xout is a cluster index image *per time*; i.e. each Xout(:,t) gives cluster labels,
% like as output from cluster --oindex option; in particular, we do not require that cluster
% indicies are unique over time on input.
% On output, Xout contains the cluster sizes (clustextent) or t-values
% (clustmass) of all clusters merged over time
%nci = max(Xin(:,1))+1; % Next cluster index  (*not* number of clusters)
for t = 2:nT
    I=Xin(:,t)>0; % mask of clusters for time t
    tclust = unique(Xin(I,t))';
    if ~isempty(tclust)
        for c = 1:length(tclust)
            idx = find(Xin(:,t)==tclust(c));   % Voxels in cluster c at time t
            Xb  = Xin(idx,t-1);        % Same voxels, but cluster indicies at time t-1
            uXb = unique(Xb(Xb>0)); % Unique labels in time t-1; assume that uXb is sorted
            
            if length(uXb)==1
                % Exactly one cluster corresponds: the time t cluster takes on t-1's index.
                % Tag it negative, to indicate this cluster in time t is newly labeled (to avoid
                % conflict between cluster indcies in times 1:t-1 and time t).
                %Xout(idx,t) = -uXb;
                idxx = logical(zeros(nV,nT));
                idxx(Xin(:,1:t-1) == uXb) = 1; % get indices of voxels from timepoint 1 to t-1 that need relabeling
                nwlab = uXb + tclust(c); % new label is the sum of all values
                Xin(idxx) = nwlab; % relabel connected clusters at previous timepoints
                Xin(idx,t) = nwlab; % relabel cluster at current timepoint
            elseif length(uXb)>1
                % Multiple clusters corresponds, so we need to merge
                idxx = logical(zeros(nV,nT));
                nwlab = sum(uXb) + tclust(c); % new label is the sum of all values
                for i=1:length(uXb)
                    idxx(Xin(:,1:t-1) == uXb(i)) = 1; % get indices of voxels from timepoint 1 to t-1 that need relabeling
                end
                Xin(idxx) = nwlab; % relabel connected clusters at previous timepoints
                Xin(idx,t) = nwlab; % relabel cluster at current timepoint
            end
        end
    end
end
toc
end


