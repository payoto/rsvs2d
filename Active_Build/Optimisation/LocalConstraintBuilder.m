function [loop]=LocalConstraintBuilder(structin,isPlot)
    
    if nargin<2
       isPlot=false; 
    end
    
    if ischar(structin)
        [structin]=PredefCases(structin);
    end
    
    loop=repmat(struct('coord',[],'isccw',true,'isinternal',false),[1,numel(structin)]);
    for ii = 1:numel(structin)
        
        switch structin(ii).name
            case 'square'
                [coord]=ProcessSquare(structin(ii));
            case 'circle'
                [coord]=ProcessCircle(structin(ii));
            case 'brick'
                [coord]=ProcessBrick(structin(ii));
                
            case 'triangle'
                [coord]=ProcessTriangle(structin(ii));
                
            otherwise
                warning('Unrecognised shape, skipping')
        end
        loop(ii).coord=coord;
    end
    if isPlot
        PlotLoop(loop,'coord')
    end
    [isInternal]=FindInternalLoop(loop,'coord');
    for ii=1:numel(loop)
        loop(ii).isinternal=isInternal(ii);
        loop(ii).isccw=CCWLoop(loop(ii).coord);
    end
end


function [structin]=PredefCases(charin)
    
    switch charin
        case '1c2b_v1'
            % design space for 1c2b cases: [0 -0.15 ;1 -0.15 ; 1 0.15 ;0 0.15]
            ii=1;
            structin(ii).name='circle';structin(ii).x=0.7;structin(ii).y=0;structin(ii).r=0.035;
            ii=2;
            structin(ii).name='brick';structin(ii).x=0.4;structin(ii).y=0.08;structin(ii).xr1=1;
            structin(ii).yr1=1;structin(ii).xr2=0.2;structin(ii).yr2=0.05;structin(ii).rot=pi/4;
            ii=3;
            structin(ii).name='brick';structin(ii).x=0.4;structin(ii).y=-0.10;structin(ii).xr1=1;
            structin(ii).yr1=1;structin(ii).xr2=0.2;structin(ii).yr2=0.05;structin(ii).rot=0;
        case '1c2b_v2'
            ii=1;
            structin(ii).name='circle';structin(ii).x=0.2;structin(ii).y=0;structin(ii).r=0.035;
            ii=2;
            structin(ii).name='brick';structin(ii).x=0.55;structin(ii).y=0.08;structin(ii).xr1=1;
            structin(ii).yr1=1;structin(ii).xr2=0.2;structin(ii).yr2=0.05;structin(ii).rot=pi/4;
            ii=3;
            structin(ii).name='brick';structin(ii).x=0.55;structin(ii).y=-0.10;structin(ii).xr1=1;
            structin(ii).yr1=1;structin(ii).xr2=0.2;structin(ii).yr2=0.05;structin(ii).rot=0;
            
        case '1c2b_v3'
            ii=1;
            structin(ii).name='circle';structin(ii).x=0.25;structin(ii).y=0;structin(ii).r=0.035;
            ii=2;
            structin(ii).name='brick';structin(ii).x=0.55;structin(ii).y=0.05;structin(ii).xr1=1;
            structin(ii).yr1=1;structin(ii).xr2=0.2;structin(ii).yr2=0.05;structin(ii).rot=pi/4;
            ii=3;
            structin(ii).name='brick';structin(ii).x=0.55;structin(ii).y=-0.06;structin(ii).xr1=1;
            structin(ii).yr1=1;structin(ii).xr2=0.2;structin(ii).yr2=0.05;structin(ii).rot=0;
    end
    
end

function [coord]=ProcessCircle(structin)
    func=@(x,y,r,pts) [(cos(linspace(0,2*pi,pts))*r+x)',(sin(linspace(0,2*pi,pts))*r+y)'];
    
    
    coord=func(structin.x,structin.y,structin.r,200);
    
end

function [coord]=ProcessSquare(structin)
    
    coord=[-1 -1; 1 -1 ;1 1; -1 1]/2*structin.l;
    
    coord=[coord(:,1)+structin.x,coord(:,2)+structin.y];
    
end

function [coord]=ProcessBrick(structin)
    
    coord=[-1 -1; 1 -1 ;1 1; -1 1]/2;
    
    coord=[coord(:,1)*structin.xr1,coord(:,2)*structin.yr1];
    d=structin.rot;
    coord=([cos(d),-sin(d);sin(d),cos(d)]*coord')';
    coord=[coord(:,1)*structin.xr2,coord(:,2)*structin.yr2];
    coord=[coord(:,1)+structin.x,coord(:,2)+structin.y];
end

function [coord]=ProcessTriangle(structin)
    
    coord=structin.coord;
end

