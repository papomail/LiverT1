%% to read .DCM volumes from the mMR Biograph
% If files are in .IMA format, either use dicom_sort_convert_main for
% conversion to DCM or radiant export.
%%
    
function [data, folderfiles, folders] = reader_Aug2018(the_folder,message)



folders=uipickfiles('FilterSpec',the_folder , 'REFilter','spirTSE|TR|0.nii\','Prompt',message,'Type', {'*.nii','NIFTI files'});

NumberLoaded=numel(folders);
clear files folderfiles
%% check if folders are already folderfiles
for ii=1:NumberLoaded
 
folderfiles{ii} = [];
[a, ~, c] = fileparts(folders{ii});
    if ~isempty(c) %If you choose file(s,) it wont look though all the directories 
        folderfiles{1} = folders;
        NumberLoaded = 1;
    end
end
%%

for ii=1:NumberLoaded 
    
        
    if isempty(folderfiles{ii})     %(NIFTI)
        folderfiles{ii}=dirrec(folders{ii},'.nii '); %looks for NIFTI
    end
    
    if isempty(folderfiles{ii}) %If still empty
       folderfiles{ii}=dirrec(folders{ii},'.img ');% look for ANALYZE
    end
    
    
     %dicom format
        if isempty(folderfiles{ii})
           folderfiles{ii}=dirrec(folders{ii},'.dcm ');
           Folder_niftii = 'NIFTI';
           [path_nii] = dcm2nii_recursive(folders,'',Folder_niftii); % convert the data to .nii
          
           message = 'Now please choose .nii files to fit. (ie. IR, mFA or LLocker)';
           folderfiles{ii}=uipickfiles('FilterSpec',[path_nii{1},Folder_niftii] , 'REFilter','.nii$','Prompt',message);
        end
        

    
    if isempty(folderfiles{ii}) %any other format (IMA)
       display('There are no files I recognise. Please if you have data as .IMA convert them to .DCM first')
    end
    
    %%Now that are in Niftii READ THEM
    
   Number_acquisitions=numel(folderfiles{ii});
for ij = 1:Number_acquisitions
    path_toread=folderfiles{ii};
    %Data_str=load_Rnii(path_toread{ij});
    
        %data{ij,ii}=load_untouch_nii(path_toread{ij});
         data{ij,ii}=nii_tool('load',path_toread{ij});

    
             
    %%%%
end

end

    
end

%% 