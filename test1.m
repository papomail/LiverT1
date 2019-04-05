function test1(varargin)



use_schedule=0; %(do not use FLIRT schedule by default)
if numel(varargin)<1
    message='Select the .nii files to do motion correction on, putting Reference image first and the Mask last. ';
    files=uipickfiles('FilterSpec',the_folder , 'REFilter','spirTSE|TR|0.nii\','Prompt',message,'Type', {'*.nii','NIFTI files'});
else
    options = struct(varargin{:});

    files={options.files};
    
    if isfield(options,'schedule')
        use_schedule=1;
        schedule=options.schedule;
    end
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% start MOCO
disp('VOLUME MOCO in progresss... please be patient');
mocotic=tic;



[mydir,~,~]=fileparts(files{1});
NumberLoaded=numel(files);

%split_ref=['cd ',dir,'; fslsplit ',refim,' aREF -z' ] ;

% Mask file
mask=files{NumberLoaded};

% image files
all_im=files(1:end-1);

%% STEP: choose one of the dimension (water) if mDixon data
clear waterfiles
splittic=tic;
for ii=1:numel(all_im)
    
    favol=files{ii};
    [~,b,~]=fileparts(favol);
    waterfiles{ii}=['W_',b,'.nii.gz'];
    disp(['Spliting ', favol])
    cmd_slice_mov=['cd ',mydir,'; fslsplit ',favol, '; rm -f vol0000.nii.gz ; rm -f vol0002.nii.gz ; rm -f vol0003.nii.gz ; mv vol0001.nii.gz ',waterfiles{ii}] ;
    %system(cmd_slice_mov);
    
end
splittime=num2str(round(toc(splittic)/60,1));
disp(['It took ',splittime ,' minutes to split the volumes.',]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Volume Motion Correction (for mFA data)
[~,r,~]=fileparts(waterfiles{1});%% reference
refim1=waterfiles{1};
for ii=1:numel(waterfiles)-1
    movim1=waterfiles{ii+1};
    [~,b,~]=fileparts(movim1);
        outim1=[b,'_3DREG'];
    moco3dtic=tic;
    disp(['Performing 3D moco of ', b, ' using ', r, ' as reference.'])
    if use_schedule
      cmd_3dmoco=['cd ',mydir,'; flirt -ref ',refim1,' -in ',movim1,' -refweight ',mask,' -inweight ',mask,' -out ',outim1, ' -schedule ',schedule] ;
    else
      cmd_3dmoco=['cd ',mydir,'; flirt -ref ',refim1,' -in ',movim1,' -refweight ',mask,' -inweight ',mask,' -out ',outim1] ;
    end
    %cmd_3dmoco=['flirt -ref ',refim1,' -in ',movim1,' -refweight ',mask1,' -inweight ',mask1,' -out ',outim1,' -cost normmi -schedule /Users/patxi/Desktop/ytransonly.sch'] ;
   % system(cmd_3dmoco);
    moco3dtime=num2str(round(toc(moco3dtic)/60,1));
    disp([b ' done in ',moco3dtime ,' minutes. (',char(datetime),').']);
    %disp(datetime);
end
