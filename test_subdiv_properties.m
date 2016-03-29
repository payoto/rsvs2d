function []=test_subdiv_properties(subMaskCell,pointsCell,eigenMat,nRefine,lvlAna)
    
    [subMaskTot]=CalcSubMaskTot(subMaskCell,nRefine,lvlAna);
    try 
        test_Basis(subMaskTot,eigenMat)
    catch ME
        warning(ME.getReport)
    end
    
    try 
        test_Frequency(subMaskTot,eigenMat)
    catch ME
        warning(ME.getReport)
    end
    
    try
        test_Projection(subMaskTot,eigenMat,pointsCell,nRefine,lvlAna)
    catch ME
        warning(ME.getReport)
    end
    
    
end


function [subMaskTot]=CalcSubMaskTot(subMaskCell,nRefine,lvlAna)
    subMaskTot=1;
    for ii=lvlAna:(nRefine)
        subMaskTot=subMaskCell{ii}*subMaskTot;
        
    end
    
end

function []=test_Basis(subMaskTot,eigenMat)
    
    figure('name','Subdivision Basis Function')
    subplot(1,2,1)
    plot(subMaskTot);
    title('Final Mask')
    
    subplot(1,2,2)
    plot(eigenMat*subMaskTot)
    title('Eigen Matrix projection')

    
end

function []=test_Frequency(subMaskTot,eigenMat)
    
    figure('name','Subdivision FFT')
    fftTarget=subMaskTot(:,3);
    subplot(1,2,1)
    plot(fftTarget);
    title('FFT target')
    
    fftres=fft(fftTarget);
    subplot(1,2,2)
    semilogy(abs(fftres))
    title('FFT result')

    
end

function []=test_Projection(subMaskTot,eigenMat,pointsCell,nRefine,lvlAna)
    
    figure('name','Subdivision surfaces')
    subplot(1,2,1)
    hold on
    axis equal
    for ii=1:nRefine
        plot(pointsCell{ii}([1:end,1],1),pointsCell{ii}([1:end,1],2))
    end
    title('All levels of subdivision')
    
    subplot(1,2,2)

    hold on
    axis equal
    startLvl=pointsCell{lvlAna};
    endLvl=pointsCell{end};
    projectLvl=eigenMat*subMaskTot*startLvl;
    %projectFact=1/1.261436458122615; % area
    projectFact=-1/sqrt(2); %bspline
    plot(startLvl([1:end,1],1),startLvl([1:end,1],2),endLvl([1:end,1],1),endLvl([1:end,1],2),projectLvl([1:end,1],1)*projectFact,projectLvl([1:end,1],2)*projectFact)
    legend('Analysis level','Final Subdivision lvl','Projection of points')
    title('analysis level and its projection')
    figure(1)
    plot(projectLvl([1:end,1],1)*projectFact,projectLvl([1:end,1],2)*projectFact)
    
end
