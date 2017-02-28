function copyfileRobust(varargin)
    
    try
        copyfile(varargin{:})
    catch ME1
        try 
            CopyFileLong(varargin{:})
        catch ME2
            
            disp(ME1.getReport)
            disp(ME2.getReport)
            throw(ME1)
        end
    end
    
    
end