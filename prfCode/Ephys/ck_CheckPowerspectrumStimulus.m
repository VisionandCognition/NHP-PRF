% NB load a stimulus mat-file from the ephys log folder

addpath(genpath('/Users/chris/Dropbox/MATLAB_NONGIT/TOOLBOX/ez'));

ri=0;
for rfs = 30:10:90
    ri=ri+1;
    RF{ri}=[];
    for xc = 0:50:800
        for yc = 0:50:800
            % pick pixel
            p=[xc yc]; %[354 354];
            ps=rfs; %60
            
            L=[];
            % check luminance fluctuations for 1 cycle
            t=0;
            for s = 1:30
                %fprintf(['Step ' num2str(s) '\n']);
                for f=1:8
                    lpix=[];
                    for x=1:ps
                        for y=1:ps
                            %lpix=[lpix ret_vid(s).img{f}(p(2)+x,p(2)+y)];
                            lpix=[lpix stimulus.img{s}(p(2)+x,p(2)+y,f)];
                        end
                    end
                    L=[L mean(lpix)];
                    %t=[t;t+1/16];
                end
            end
            RF{ri} = [RF{ri}; L];
        end
    end
end
%%

% upsample to 100 Hz
t=0:1/15:240/15; t(end)=[];

L2=[];
t2=0:0.01:max(t);
for i=1:length(t2)
    idx=find(t<=t2(i),1,'last');
    L2=[L2;L(idx)];
end


%plot(t2,L2)

ez_powerspec(L2,100)