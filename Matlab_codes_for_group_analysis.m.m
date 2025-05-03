    clear
    warning('off')
    % Set global variables and parameters
    derivatives_dir='/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives';
    % cd(derivatives_dir)
    subjects=readtable('/home/koba/Desktop/Stroke/scripts/subjects.csv');
    behavioral=readtable('/home/koba/Desktop/Stroke/scripts/behavioral.csv');
    parcellation_nii=load_nii('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/source/yeo400_resampled.nii.gz');
    parcellation=reshape(parcellation_nii.img,[],1);
    parcellation_mask_whole=parcellation>0;
    parcellation_mask_left=parcellation>0 & parcellation <201; % first half is left hemisphere
    parcellation_mask_right=parcellation>200;
    parcellation_reduced=parcellation(parcellation_mask_whole>0);

    eye_mask_nii=load_nii('/home/koba/Desktop/brainhack/nau_mask_resampled.nii.gz');
    eye_mask=logical(reshape(eye_mask_nii.img,[],1));


    [surf_lh, surf_rh] = load_conte69();
    labeling = load_parcellation('schaefer',400);
    addpath(genpath('/media/koba/MULTIBOOT/blindness_gradients/source/toolboxes/BrainStat-master'))
    [network_names, colormap] = fetch_yeo_networks_metadata(7);
    schaefer_400 = fetch_parcellation('fsaverage5', 'schaefer', 400);
    rmpath(genpath('/media/koba/MULTIBOOT/blindness_gradients/source/toolboxes/BrainStat-master'))
    addpath(genpath("/home/koba/Desktop/kramer"))
    
    TR = 2.5;
    samplingRate = 1 / TR;
    lowFreq = 0.01; % Lower cutoff frequency in Hz
    highFreq = 0.1; % Upper cutoff frequency in Hz
    filterOrder = 5; % Filter order
    [b, a] = butter(filterOrder, [lowFreq, highFreq] * 2 / samplingRate);
    
    lags=-3:3; % TR, not seconds
    num_control=sum(subjects.subj_type);
    num_stroke=size(subjects,1)-num_control;
    
    % gather the necessary files
    func_list=char(table2array(readtable('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/source/single_runs.txt',Delimiter=',',ReadVariableNames=false)));
    mask_list=char(table2array(readtable('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/source/mask_files.txt',Delimiter=',',ReadVariableNames=false)));
    confounds_list=char(table2array(readtable('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/source/confounds.txt',Delimiter=',',ReadVariableNames=false)));
    seven_network_positions=table2array(readtable('/media/koba/MULTIBOOT/net/ascratch/people/plgkoba/stroke_BIDS_firs_sessions/derivatives/source/7netpositions.txt',ReadVariableNames=false));
    network_names={'Visual','Somatomotor','Dorsal Att.','Ventral Att.','Limbic','Control','Default Mode'};
   
%% Postprocessing: confound regression, band-pass filtering, lag calculation, saving niftis, lag correction, creating correlation matrices

for i=1:size(func_list,1)
    disp(i)
    % load the files

    func_file=load_nii(strrep(func_list(i,:),' ', ''));
    confound_file=readtable(strrep(confounds_list(i,:),' ', ''),"FileType","text");
    %     mask_file=load_nii(strrep(mask_list(i,:),' ', ''));

    % mask the functional data and keep the confounds of interest
    func_data=reshape(func_file.img,[],size(func_file.img,4));
    %mask_data=reshape(mask_file.img,[],1);
    func_masked=func_data(eye_mask,:); % func data is converted to 2D and masked
    confounds_of_interest = {'framewise_displacement' 'trans_x' 'trans_y' 'trans_z' 'rot_x' 'rot_y' 'rot_z'};

    confounds_data=[table2array(confound_file(:,confounds_of_interest)) (1:size(func_data,2))']; % adding linear trend
    confounds_data(isnan(confounds_data))=0; % change nan with zero

    % Regress te confounds our and apply band-pass filter
    func_clean=zeros(size(func_masked));
    for j=1:size(func_masked,1)
        [~,~,residuals] = regress(func_masked(j,:)',[ones(size(func_masked,2),1) confounds_data]);
        func_clean(j,:) = filter(b, a, residuals);
    end
    func_clean_norm = (func_clean - min(func_clean, [], 2)) ./ (max(func_clean, [], 2) - min(func_clean, [], 2));
     sub_id=func_list(i,10:19);
    % if strfind(func_list(i,:),'_space')-(strfind(func_list(i,:),'n-')+2)==1
    %     run=strfind(func_list(i,:),'n-')+2;
    % else
    %     run=[strfind(func_list(i,:),'n-')+2 strfind(func_list(i,:),'n-')+3];
    % end
    run=0;
    % Save the image
    template=func_file;
    template_data=zeros(size(func_data));
    template_data(eye_mask,:)=func_clean_norm;
    template.img=reshape(template_data,[size(func_file.img)]);
    
    file_name='brainhack/preprocessed_files/only_mc/SUBID_run-RUN_eyes_only_mc.nii.gz';
    file_name=strrep(file_name,'SUBID',sub_id);
    file_name=strrep(file_name,'RUN',string(run));
    save_nii(template,char(file_name))


    % % % Smooth the image - gunzip it first
    % file_path_gunzip='brainhack/preprocessed_files/temp';
    % gunzip(file_name,file_path_gunzip)
    % 
    % file_name_gunzip='brainhack/preprocessed_files/temp/SUBID_run-RUN_eyes_notsmoothed.nii';
    % file_name_gunzip=strrep(file_name_gunzip,'SUBID',sub_id);
    % file_name_gunzip=strrep(file_name_gunzip,'RUN',func_list(i,run));
    % 
    % 
    % 
    % matlabbatch{1}.spm.spatial.smooth.data = arrayfun(@(x) sprintf('%s,%d', file_name_gunzip, x), 1:128, 'UniformOutput', false)';
    % matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
    % matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    % matlabbatch{1}.spm.spatial.smooth.im = 0;
    % matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    % spm('defaults', 'fmri'); % Load default settings for fMRI
    % spm_jobman('initcfg'); % Initialize the SPM job manager
    % 
    % spm_jobman('run', matlabbatch)

   
end

eye_mask_3d=eye_mask_nii.img;
eye_mask_right=eye_mask_3d;
eye_mask_right(1:24,:,:)=0;
eye_mask_left=eye_mask_3d;
eye_mask_left(24:end,:,:)=0;

template=eye_mask_nii;
template.img=eye_mask_right;
save_nii(template,char('brainhack/nau_mask_right.nii.gz'))

template=eye_mask_nii;
template.img=eye_mask_left;
save_nii(template,char('brainhack/nau_mask_leftt.nii.gz'))



eye_mask_inverted=load_nii('eye_mask_zosia.nii');
eye_mask_inverted_data=eye_mask_inverted.img;

eye_mask_inverted.img=abs(eye_mask_inverted_data-1);
save_nii(eye_mask_inverted,char('eye_mask_zosia_fixed.nii'))
