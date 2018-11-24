function [corrp tstats] = osl_clustertf(c,thresh,nP,con,varcope_time_smooth_std,tres,permmeth,varargin)

% [corrp tstats] = osl_clustertf(c,thresh,nP,con,varcope_time_smooth_std,tres)
%
% Only works for a group average at the moment to do permutation tests 
% (using sign flipping).
%
% c is data where:
% nS = size(c,1); %number of 'subjects' (recording sites)
% nF = size(c,2); %number of frequencies
% nT = size(c,3); %number of timebins
%
% thresh is threshold to use on c
%
% nP is num of permutations
%
% con is the connectivity criterion. Could be 6(surface)
%               18(edge) or 26(corner). For the 2D case here, these will
%               correspond to 4, 8 and 8 respectively.
%
% tres is temporal res of data in secs
%
% varcope_time_smooth_std is smoothing to do on time courses of group
% varcopes
%
% LH and MWW 2012

if nargin<2
  thresh = 3.1;
end

if nargin<3
  nP = 500;
end

if nargin<4
  con = 26;
end

if nargin<5
  varcope_time_smooth_std = 0;
  tres=[];
end


nS = size(c,1); %number of 'subjects' (recording sites)
nF = size(c,2); %number of frequencies
nT = size(c,3); %number of timebins

if length(size(c))==4
    c = squeeze(c(:,:,:,varargin{1,1}));
end

cr = reshape(c,nS,nF*nT);

% setup gauss for temporal smoothing  
if(varcope_time_smooth_std>0)
    fft_gauss=fftshift(gauss(varcope_time_smooth_std/tres,1,nT)');
end;

%build up a null distribution of the maximum cluster size for each
%permuatation i:
nulldist=zeros(nP,1);

for i = 1:nP

    dm = ((rand(nS,1)>0.5)-0.5)*2;
  
    [cg,vg,tg] = ols(cr,dm,1);

    % variance smooth over time:
    if(varcope_time_smooth_std>0)
        vg=reshape(vg,nF,nT);

        for f=1:size(vg,1),
            dat=permute(vg(f,:),[2, 1, 3, 4]);
            dat2 = fftconv(dat,fft_gauss);
            vg(f,:)=dat2;
        end;
        
        cg=reshape(cg,nF,nT);
        tg=cg./sqrt(vg);
        
    else
        tg = reshape(tg,nF,nT);
    end;
    
    for jj=1:2
        if jj == 1
            [imlabel LL] = spm_bwlabel(double(tg>thresh),con); %positive clusters
        elseif jj == 2
            [imlabel LL] = spm_bwlabel(double(tg<-thresh),con); %negative clusters
        end
        tmp = unique(imlabel);
        tmp(tmp==0) = [];
        nL = 0;
        
        if ~isempty(tmp)
            for kk = 1:numel(tmp); %loop over clusters
                k = tmp(kk);
                if strcmp(permmeth,'clustextent')
                    nL(k) = sum(sum(imlabel==k));
                elseif strcmp(permmeth,'clustmass')
                    nL(k) = sum(tg(imlabel==k));
                end
                tsmax(jj)=max(nL); %two-sided max
            end
            nulldist(i)=max(abs(tsmax)); 
        end
    end
end

%run a one sample t-test on data and compare cluster size to
%permutation data
dm = ones(nS,1);
[cg,vg,tg] = ols(cr,dm,1);

% variance smooth over time:
if(varcope_time_smooth_std>0)
    vg=reshape(vg,nF,nT);

    for f=1:size(vg,1),
        dat=permute(vg(f,:),[2, 1, 3, 4]);
        dat2 = fftconv(dat,fft_gauss);  
        vg(f,:)=dat2;
    end;

    cg=reshape(cg,nF,nT);
    tg=cg./sqrt(vg);

else
    tg = reshape(tg,nF,nT);        
end;

for ll=1:2
    if ll == 1
        [imlabel LL] = spm_bwlabel(double(tg>thresh),con); %positive clusters
        cpimlabel = zeros(size(imlabel));
    elseif jj == 2
        [imlabel LL] = spm_bwlabel(double(tg<-thresh),con); %negative clusters
    end
    tmp = unique(imlabel);
    sprintf('%d clusters found',size(tmp,1));
    tmp(tmp==0) = [];
    if ~isempty(tmp)
        for kk = 1:numel(tmp) %loop over clusters
            k = tmp(kk);
            if strcmp(permmeth,'clustextent')
                nL = abs(sum(sum(imlabel==k)));
            elseif strcmp(permmeth,'clustmass')
                nL= abs(sum(tg(imlabel==k)));
            end
            cp = mean(nL>nulldist);
            if cp == 1 % if no cluster is larger
                cp = 1-1/nP; % set p-value to 1/number of permutations
            end
            cpimlabel(imlabel==k)=cp;
            sprintf('Cluster %d: value %.2f,  p = %.4f',kk,double(nL),(1-cp));
        end
    end
end
corrp = cpimlabel;
tstats=tg;
