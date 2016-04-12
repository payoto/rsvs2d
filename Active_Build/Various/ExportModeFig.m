function ExportModeFig(figStart,nFig,casestr)
    folderName=['..\results\Development Notes\Profile Modes\',casestr];
    mkdir(folderName);
    for ii=1:nFig
        fileName=[folderName,'\mode',...
            num2str(ii)];
        print(ii+figStart-1,[fileName,'.png'],'-dpng','-r200');
        
    end

    
end