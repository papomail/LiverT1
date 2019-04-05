EXT1files=FilterFiles('/Users/patxi/Desktop/REVITA22','EXTREMEsmoothB1.mat');
NN=numel(EXT1files);
clear ScanID FolderNum
for ii=1:NN
    newfolder=['/Users/patxi/Desktop/REVITA2_T1MAPS/',num2str(NN+1-ii)];
    mkdir(newfolder);
    matfile=EXT1files(ii).name;
    copyfile(matfile,newfolder)
    
    niftifile=[matfile(1:end-3),'nii.gz'];
        copyfile(niftifile,newfolder)
        
        FolderNum{ii}=num2str(NN+1-ii);
        ScanID{ii}=fileparts(fileparts(matfile));
end


T =table(FolderNum(:),ScanID(:),'VariableNames',[{'FolderNum'},{'ScanID'}]);

writetable(T,'/Users/patxi/Desktop/REVITA2_T1MAPS/NamesCode.xls')


