function [ dinfo ] = dparse_folders(din, filter)
% *Modified version of dparse to read files inside the input folders that match the given name-filer.
% *Assign RealWorldValueSlope=1 and RealWorldValueIntercept=0 if RealWorldValueMappingSequence exists but is empty
% *Use data_process_noASL instead of data_process

% Patxi Torrealdea



%DPARSE DICOM file parse (SingleFrame DICOM)
%   Parses single-frame DICOM files within a directory structure. Calls
%   cast_private to fix some PACs and ransfer issues
% dinfo = dparse
%
% dinfo is returned with key fields, including an integer slice index in
% field sl
%
% See also D2MAT XMLPARSE CAST_PRIVATE DATA_PROCESS
%
% David Atkinson
% $Id: dparse.m 341 2012-06-17 21:23:24Z ucacdat $

global VOLPATH   % ensures paths persist across calls to this 
                 % function (for ease of use and testing, not essential)
  

    VOLPATH = din ;


% These DICOM 3 tags will be placed in dinfo if present in files
tags = {'Filename', 'PixelSpacing', 'ImageOrientationPatient', ...
    'ImagePositionPatient', 'Height', 'Width', 'SeriesNumber', ...
    'SliceThickness', 'AcquisitionTime', 'SeriesTime', ...
    'TemporalPositionIdentifier', ...       % Philips (helpful)
    'ImageType', ...  % Useful for Dixon
    'InversionTime',...
    'TriggerTime', ...
    'FlipAngle', 'Private_2001_1023', ... % 1st can lose precision
    'RepetitionTime', 'Private_2005_1030', ...  % 1st can be zero if small!
    'DiffusionBValue', 'DiffusionGradientOrientation', ...
    'Private_0019_100c', ...  % Siemens b-value
    'ProtocolName', ...
    'FrameOfReferenceUID', ...
    'MRAcquisitionType',...
    'AcquisitionNumber',...
    'EchoNumber', ...
    'EchoTime',...
    'RescaleSlope','RescaleIntercept', 'Private_2005_100e', ...
    'RealWorldValueMappingSequence', ...
    'Private_2005_1429',... % MJS adding the label type for single frame images
    'Private_2001_1008',... % MJS adding the phase number for multiphase images
    'Private_2001_1017',... % MJS adding the number of phases for multiphase images
    'Private_2005_1409','Private_2005_140a', ... % DA RescaleSlope and Intercept
} ;


%disp('Select folder with DICOM files.')
%voldir = uigetdir(VOLPATH,'Select folder containing single-frame DICOM files.') ;

voldir=VOLPATH

if voldir == 0 ; return; end ;
VOLPATH = voldir ;

% Assemble d, a list of directories. Include sub-folders one level down.


disp(voldir);
dtoptop = dir(voldir) ;




realdtop={};
counter=0;
for ii=1:numel(dtoptop); 
    if dtoptop(ii).isdir && ~isempty(strfind(dtoptop(ii).name,filter))
        counter=counter+1;
        realdtop{counter}=dtoptop(ii);
    end
end
dtop=realdtop;

if numel(dtop)==0
    disp(['Cannot find any folder matching the filter. Please review folder ',voldir, ' and check the names match the filter given.'])
    dinfo=0;
    return 
end

isubf = 0 ;
d = dtop ;

for ident = 1:numel(d)
    if dtop(ident).isdir 
        if strcmp(dtop(ident).name,'.')==0 && ...
                strcmp(dtop(ident).name,'..')==0
            isubf = isubf + 1;
            subf{isubf} = dtop(ident).name ;
            %Produce directory structure for each subfolder, 
            % add sub-folder to name and append overall direc 
            dsubf = dir(fullfile(voldir,subf{isubf})) ;
            for ident = 1:numel(dsubf)
                dsubf(ident).name = fullfile(subf{isubf}, dsubf(ident).name) ;
            end
            d = [ d ; dsubf ];
        end
    end
end

ikeep = [] ;

% Select dir entries that are DICOM files
% Later exclude the PS_* and XX_* files output by ViewForum, DICOMDIR
% and other DICOM compatible files that are not underlying acquired images.
% (These require the info structure to be read)

for ident = 1:numel(d)
    if ~d(ident).isdir && d(ident).bytes > 128
        if isdicom(fullfile(voldir,d(ident).name))
            ikeep = [ikeep ident] ;
        end
    end
end

d = d(ikeep) ;
nd = numel(d) ;

if exist('subf','var')
    nsubf = numel(subf) ;
else
    nsubf = 0 ;
end

df = 1 ; % field counter

% dinfo = struct('RescaleSlope',num2cell(ones([1 nd])), 'RescaleIntercept', ...
%     num2cell(zeros([1 nd])), 'Private_2005_100e',num2cell(ones([1 nd]))) ;
dinfo = struct ;

npsxx = 0 ;

hw = waitbar(0,['Processing ',num2str(nd),' files in ',...
    num2str(nsubf),' folder(s)']) ;

for id = 1:nd
    waitbar(id/nd,hw) ;
    info = dicominfo(fullfile(voldir,d(id).name)) ;
    
    % insert skip if View Forum XX or PS files, DICOMDIR, non-geom
    if ~isfield(info,'DirectoryRecordSequence') && ...  % DICOMDIR
            isfield(info,'ImagePositionPatient') && ...  % 
            strcmp(info.SOPClassUID, '1.2.840.10008.5.1.4.1.1.11.1') == 0 && ...
            strcmp(info.SOPClassUID, '1.2.840.10008.5.1.4.1.1.66') == 0
         
             % 1.2.840.10008.5.1.4.1.1.11.1 Grayscale Softcopy Presentation State
             % 1.2.840.10008.5.1.4.1.1.66   Raw data
        for itag = 1:length(tags)
            if isfield(info,tags{itag})
                if strcmp(tags{itag},'AcquisitionTime') || strcmp(tags{itag},'SeriesTime')
                    dinfo = setfield(dinfo,{df},tags{itag},str2num(getfield(info,tags{itag}))) ;
                else
                    dinfo = setfield(dinfo,{df},tags{itag},getfield(info,tags{itag})) ;
                end
            else
                % warning([tags{itag},' is not in info'])
            end
        end
%       if ~isfield(info,'Private_2005_100e')
%           % Not in input, default value of 1 is best replaced with 1/RescaleSlope
%           dinfo(df).Private_2005_100e = 1/dinfo(df).RescaleSlope ;
%       end
        if isfield(dinfo(df), 'RealWorldValueMappingSequence')
            if ~isempty(dinfo(df).RealWorldValueMappingSequence)
            dinfo(df).RealWorldValueIntercept = dinfo(df).RealWorldValueMappingSequence.Item_1.RealWorldValueIntercept;
            dinfo(df).RealWorldValueSlope = dinfo(df).RealWorldValueMappingSequence.Item_1.RealWorldValueSlope;
            else
              dinfo(df).RealWorldValueIntercept = [0];
              dinfo(df).RealWorldValueSlope =[1];
            end
        end

        df = df + 1 ;
          
        
    else % XX or PS file test OR DICOMDIR
        npsxx = npsxx + 1 ;
    end
end % id loop

close(hw)  % close waitbar

ndf = df - 1 ; 
dinfo(ndf+1:end) = [] ;
if npsxx > 0 
    str = [' (ignored ',num2str(npsxx),' files).'];
else
    str = ['.'];
end
disp(['Processed ',num2str(ndf),' DICOM files', str])

if isfield(dinfo,'RealWorldValueMappingSequence')
    % if ~isempty(dinfo(df).RealWorldValueMappingSequence)
     dinfo = rmfield(dinfo,'RealWorldValueMappingSequence') ;
    % end
end

dinfo = cast_private(dinfo) ; % Fixes unknown value representaton 
                              % for Private fields in some circumstances
%dinfo = data_process(dinfo) ;
dinfo = data_process_noASL(dinfo) ;
end


