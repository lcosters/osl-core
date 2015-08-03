function glean_model(GLEAN)
% Runs the model stage of GLEAN

if ~exist(GLEAN.model,'file')
    
    % Concatenate data:
    dataConcat = cell(1,numel(GLEAN.data)); % cell arrays grow better than arrays
    subIndx    = cell(1,numel(GLEAN.data)); % cell arrays grow better than arrays
    for session = 1:numel(GLEAN.data)
        D = spm_eeg_load(GLEAN.data(session).subspace);
        %if isfield(GLEAN.settings.envelope,'freqbands') && numel(GLEAN.settings.envelope.freqbands) > 1
        %    % rearrange channels x frequencies as [c1f1, ... ,c1fn, c2f1, ...]
            dataConcat{session} = reshape(D(:,:,:,:),[],D.nsamples,D.ntrials);
        %else
        %    dataConcat{session} = D(:,:,:);
        %end
        
        % Commented code above should be deleted once 4D MEEG object is
        % used consistently
        
        subIndx{session} = session*ones(1,D.nsamples);

    end
    dataConcat = cell2mat(dataConcat);
    subIndx = cell2mat(subIndx); %#ok
    
    dataConcat = normalise(dataConcat,2); % MOVE THIS SOMEWHERE BETTER!
    
    switch lower(char(fieldnames(GLEAN.settings.model)))
        case 'hmm'
            hmm = osl_hmm_infer(dataConcat,struct('K',GLEAN.settings.model.hmm.nstates,'order',0,'Ninits',GLEAN.settings.model.hmm.nreps,'zeromean',0)); %#ok
            save(GLEAN.model,'hmm','subIndx')
        case 'ica'
            nICs = GLEAN.settings.model.ica.order;
            [ica.tICs,ica.SM,~] = fastica(dataConcat,'g','tanh','lastEig',nICs,'numOfIC',nICs,'approach','symm');
            save(GLEAN.model,'ica','subIndx')
    end
    
end

end