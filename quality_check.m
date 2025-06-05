% Load the quality control file and prepare the file paths
clear
warning('off')
qa_table=readtable('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS/stroke_dataset_quality_control.csv');
bids_root={'/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/'};
anat_root={'stroke_BIDS_NUMBER_sessions/derivatives/fmriprep/sub-SUBID/ses-SESSION/anat/sub-SUBID_ses-SESSION_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii'};
func_root={'stroke_BIDS_NUMBER_sessions/derivatives/fmriprep/sub-SUBID/ses-SESSION/func/sub-SUBID_ses-SESSION_task-rest_run-RUNORDER_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz'};
confounds_root={'stroke_BIDS_NUMBER_sessions/derivatives/fmriprep/sub-SUBID/ses-SESSION/func/sub-SUBID_ses-SESSION_task-rest_run-RUNORDER_desc-confounds_timeseries.tsv'};

% Load the toolboxes
addpath(genpath('/media/koba/MULTIBOOT/blindness_gradients/source/toolboxes/nifti'))

% Get the confounds file
confounds_of_interest = {'framewise_displacement'  'csf' 'white_matter' 'trans_x' 'trans_y' 'trans_z' 'rot_x' 'rot_y' 'rot_z' ...
    'csf_derivative1' 'white_matter_derivative1' 'trans_x_derivative1' 'trans_y_derivative1' 'trans_z_derivative1' 'rot_x_derivative1' 'rot_y_derivative1' 'rot_z_derivative1' ...
    'csf_derivative1_power2' 'white_matter_derivative1_power2' 'trans_x_derivative1_power2' 'trans_y_derivative1_power2' 'trans_z_derivative1_power2' 'rot_x_derivative1_power2' 'rot_y_derivative1_power2' 'rot_z_derivative1_power2' ...
    'csf_power2' 'white_matter_power2' 'trans_x_power2' 'trans_y_power2' 'trans_z_power2' 'rot_x_power2' 'rot_y_power2' 'rot_z_power2'};

% Parcellation mask
parcellation_nii=load_nii('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_first_sessions/derivatives/source/yeo400_resampled.nii.gz');
parcellation=reshape(parcellation_nii.img,[],1);
parcellation_mask_whole=parcellation>0;
parcellation_mask_left=parcellation>0 & parcellation <201; % first half is left hemisphere
parcellation_mask_right=parcellation>200;
parcellation_reduced=parcellation(parcellation_mask_whole>0);

% Parameters for band-pass filtering
TR = 2;
samplingRate = 1 / TR;
lowFreq = 0.01; % Lower cutoff frequency in Hz
highFreq = 0.1; % Upper cutoff frequency in Hz
filterOrder = 5; % Filter order
[b, a] = butter(filterOrder, [lowFreq, highFreq] * 2 / samplingRate);

lags=-3:3; % TR, not seconds

% Load the metdata
participants=readtable('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS/participants.tsv','FileType','text','Delimiter', '\t');

% Get the reference gradient

[surf_lh, surf_rh] = load_conte69();
labeling = load_parcellation('schaefer',400);
addpath(genpath('/media/koba/MULTIBOOT/blindness_gradients/source/toolboxes/BrainStat-master'))
[network_names, colormap] = fetch_yeo_networks_metadata(7);
schaefer_400 = fetch_parcellation('fsaverage5', 'schaefer', 400);
rmpath(genpath('/media/koba/MULTIBOOT/blindness_gradients/source/toolboxes/BrainStat-master'))
addpath(genpath("/home/koba/Desktop/kramer"))

conn_matrices = load_group_fc('schaefer',400);
reference_gradient = GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);
m = (conn_matrices.schaefer_400);
m=atanh(m);
m(isinf(m))=1;
reference_gradient = reference_gradient.fit(m);
gradient_func=GradientMaps('kernel','na','approach','dm','alignment','pa','n_components',10);


% Initiate the structure
stroke_structure(size(qa_table,1),1)=struct();
for i=1:size(qa_table,1)
    
    % Add the first info to the structure
    temp_struct = table2struct(qa_table(i,:), 'ToScalar', true);
    field_names = fieldnames(temp_struct);
    % Assign fields one by one into stroke_structure(i)
    for f = 1:length(field_names)
        stroke_structure(i).(field_names{f}) = temp_struct.(field_names{f});
    end

    % Create the specific paths for the files
    SUBID=qa_table(i,:).ID;
    SESSION=qa_table(i,:).Session;
    run=qa_table(i,:).Run;

    if strcmp(SESSION,'control') || strcmp(SESSION,'acute')
        NUMBER={'first'};
    elseif strcmp(SESSION,'control2') || strcmp(SESSION,'followup')
        NUMBER={'second'};
    elseif strcmp(SESSION,'followup2')
        NUMBER={'third'};
    end
    disp([SUBID SESSION])


    anat_file=strcat(bids_root,anat_root);
    anat_file=strrep(anat_file,'SESSION',SESSION);
    anat_file=strrep(anat_file,'SUBID',SUBID);
    anat_file=strrep(anat_file,'NUMBER',NUMBER);

    func_file=strcat(bids_root,func_root);
    func_file=strrep(func_file,'SESSION',SESSION);
    func_file=strrep(func_file,'SUBID',SUBID);
    func_file=strrep(func_file,'NUMBER',NUMBER);

    confound_file=strcat(bids_root,confounds_root);
    confound_file=strrep(confound_file,'SESSION',SESSION);
    confound_file=strrep(confound_file,'SUBID',SUBID);
    confound_file=strrep(confound_file,'NUMBER',NUMBER);
    
    run_array=str2double(strsplit(run{1}, '-'));

    corrmats=zeros(max(parcellation_reduced),max(parcellation_reduced),size(run_array,2));
    corrmats_corrected=zeros(max(parcellation_reduced),max(parcellation_reduced),size(run_array,2));

    gradients_corrmat=zeros(max(parcellation),10,size(run_array,2));
    lambda_corrmat=zeros(19,size(run_array,2));
    gradients_corrmat_corrected=zeros(max(parcellation),10,size(run_array,2));
    lambda_corrmat_corrected=zeros(19,size(run_array,2));


    % Go over the functional files
    for j=1:size(run_array,2)

        % Get the functional data and reshape it to 2D
        func_run=strrep(func_file,'RUNORDER',num2str(run_array(j)));
        func_run_data =  load_nii(func_run{1});
        func_data=reshape(func_run_data.img,[],size(func_run_data.img,4));

        % Get the confounds data
        confound_file_run=strrep(confound_file,'RUNORDER',num2str(run_array(j)));
        confounds_table=readtable(confound_file_run{1},"FileType","text");
        confounds_data=[table2array(confounds_table(:,confounds_of_interest)) (1:size(func_data,2))']; % adding linear trend
        confounds_data(isnan(confounds_data))=0; % change nan with zero

        % Regress te confounds our and apply band-pass filter
        func_masked=func_data(logical(parcellation_mask_whole),:); % func data is converted to 2D and masked
        func_clean=zeros(size(func_masked));
        for k=1:size(func_masked,1)
            [~,~,residuals] = regress(func_masked(k,:)',[ones(size(func_masked,2),1) confounds_data]);
            func_clean(k,:) = filter(b, a, residuals);
        end

        % Get the global signal -- not using fmriprep outpits because we
        % use mean signal from signal hemispheres of stroke patients

        rows_id_match = contains(participants.participant_id, SUBID);
        rows_basic_event = contains(participants.redcap_event_name, 'basic', 'IgnoreCase', true);
        final_row_basic = rows_id_match & rows_basic_event;
        basic_info = participants(final_row_basic,:);


        if strcmp(SESSION,'control')
            arm='visit_1_arm_2';
        elseif strcmp(SESSION,'control2')
            arm='visit_2_arm_2';
        elseif strcmp(SESSION,'acute')
            arm='acute_arm_1';
        elseif strcmp(SESSION,'followup')
            arm='3month_arm_1';
        elseif strcmp(SESSION,'followup2')
            arm='1year_arm_1';
        end

        rows_session = contains(participants.redcap_event_name, arm, 'IgnoreCase', true);
        final_row_session = rows_id_match & rows_session;
        session_info = participants(final_row_session,:);

        if contains(SUBID,'CON')
            parcellation_in_brainmask=parcellation_mask_whole;
            gs=mean(func_data(parcellation_in_brainmask,:));
        elseif contains(SUBID,'PAT')

            if basic_info.lesion_side == 0 % 0 is for left
                parcellation_in_brainmask=parcellation_mask_right;
                gs=mean(func_data(parcellation_in_brainmask,:));
            elseif basic_info.lesion_side == 1
                parcellation_in_brainmask=parcellation_mask_left;
                gs=mean(func_data(parcellation_in_brainmask,:));
            end
        end

        % Calculate the lag
        corrs_lag=zeros(length(lags),size(func_clean,1));
        for k=1:length(lags)
            gsr_lag=circshift(gs,lags(k));
            corrs_lag(k,:)=corr(gsr_lag',func_clean');
        end

        corrs_lag(isnan(corrs_lag))=0;
        lag_pos=zeros(1,size(func_clean,1));
        for k=1:length(corrs_lag)
            vec=(corrs_lag(:,k));
            laginfo=(find(vec==max(vec)));
            if isempty(laginfo)
                lag_pos(k)=0;
            else
                lag_pos(k)=lags(laginfo(1));
            end
        end

        % lags_all(i,:)=lag_pos;
        func_corrected=zeros(size(func_clean));
        for k=1:size(func_clean,1)
            func_corrected(k,:)=circshift(func_clean(k,:),lag_pos(k)*-1);
        end

        % Get the mean time series from each ROI
        mean_ts=zeros(max(parcellation_reduced),size(func_corrected,2));
        mean_ts_corrected=zeros(max(parcellation_reduced),size(func_corrected,2));

        for k=1:max(parcellation_reduced)
            m=func_clean(parcellation_reduced==k,:);
            mean_ts(k,:)=mean(m(any(m,2),:));

            m=func_corrected(parcellation_reduced==k,:);
            mean_ts_corrected(k,:)=mean(m(any(m,2),:));
        end

        % Get the Z transformed correlation matrices
        corrmat=atanh(corr(mean_ts'));
        corrmat(isinf(corrmat))=1;
        corrmats(:,:,j)=corrmat;

        corrmat_corrected=atanh(corr(mean_ts_corrected'));
        corrmat_corrected(isinf(corrmat_corrected))=1;
        corrmats_corrected(:,:,j)=corrmat_corrected;

        % Calculate the gradients
   

        gradient_corrmat=gradient_func.fit(corrmat,'reference',reference_gradient.aligned{1},'sparsity',90);
        gradients_corrmat(:,:,j)=gradient_corrmat.aligned{1}(:,:);
        lambda_corrmat(:,j)=gradient_corrmat.lambda{1};


        
        gradient_corrmat_corrected=gradient_func.fit(corrmat_corrected,'reference',reference_gradient.gradients{1},'sparsity',90);
        gradients_corrmat_corrected(:,:,j)=gradient_corrmat_corrected.aligned{1};
        lambda_corrmat_corrected(:,j)=gradient_corrmat_corrected.lambda{1};
    end

    % Take the mean of the correlation matrices and their gradients
    mean_corrmats=mean(corrmats,3);
    gradient_mean_corrmats=gradient_func.fit(mean_corrmats,'reference',reference_gradient.gradients{1},'sparsity',90);
    

    mean_corrmats_corrected=mean(corrmats_corrected,3);
    gradient_mean_corrmats_corrected=gradient_func.fit(mean_corrmats_corrected,'reference',reference_gradient.gradients{1},'sparsity',90);


    % Merge the data in the structure
    stroke_structure(i).corrmats=corrmats;
    stroke_structure(i).gradient_corrmat=gradients_corrmat;
    stroke_structure(i).lambda_corrmats=lambda_corrmat;

    stroke_structure(i).mean_corrmats=mean_corrmats;
    stroke_structure(i).gradient_corrmat_mean=gradient_mean_corrmats.aligned{1};
    stroke_structure(i).lambda_corrmats_mean=gradient_mean_corrmats.lambda{1};

    stroke_structure(i).corrmats_corrected=corrmats_corrected;
    stroke_structure(i).gradient_corrmat_corrected=gradients_corrmat_corrected;
    stroke_structure(i).lambda_corrmats_corrected=gradient_func;

    stroke_structure(i).mean_corrmats_corrected=mean_corrmats_corrected;
    stroke_structure(i).gradient_corrmat_corrected_mean=gradient_mean_corrmats_corrected.aligned{1};
    stroke_structure(i).lambda_corrmats_corrected_mean=gradient_mean_corrmats_corrected.lambda{1};

    temp_struct = table2struct(basic_info(:, 3:302), 'ToScalar', true);
    field_names = fieldnames(temp_struct);
    for f = 1:length(field_names)
        stroke_structure(i).(field_names{f}) = basic_info.(field_names{f});
    end


    temp_struct = table2struct(session_info(:, 303:end), 'ToScalar', true);
    field_names = fieldnames(temp_struct);
    for f = 1:length(field_names)
        stroke_structure(i).(field_names{f}) = session_info.(field_names{f});
    end
    
    whos('stroke_structure').bytes/1000000

    if mod(i, 10) == 0
        filename = sprintf('stroke_structure_iter_%d.mat', i);
        save(filename, 'stroke_structure', '-v7.3');
        fprintf('Saved iteration %d to %s\n', i, filename);
    end


end
