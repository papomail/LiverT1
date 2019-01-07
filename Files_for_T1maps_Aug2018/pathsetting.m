%% set original path
function [the_folder] = pathsetting(varargin)

   
    the_folder = pwd;
    out=regexp(the_folder,slsh,'split');

    options = varargin;


        if ismac
            the_folder = [slsh,out{2},slsh,out{3},slsh];
        elseif ispc
            the_folder = [out{1},slsh,out{2},slsh, out{3},slsh];
        end  

    if ~isempty(options)
        Newdirectory = options{1};
        the_folder = [the_folder,Newdirectory];
    end
    
end