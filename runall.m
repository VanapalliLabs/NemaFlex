
clear; clc
tic;
list = dir('*.avi*');  %save info structure for all AVI or AVI files
Nfiles = length(list);  %number of AVI to process

for l=1:Nfiles
    try
        FileAVI = list(l).name;    
        r = skel2v4n(FileAVI, 'magnification',4,'pixelsize',6.50,'diameter',40,'spacing',80,'verb',1);
        M = pl2mv(r,'thick');
        D = pl2data(r,'dist');
        %Make directory for the respective videos and move the R, M & D
        %files along with the video
        if mkdir(FileAVI(1:end-4)) == 1
            save([FileAVI(1:end-4),'/','R_',FileAVI(1:end-4)],'r')
            save([FileAVI(1:end-4),'/','M_',FileAVI(1:end-4)],'M')
            save([FileAVI(1:end-4),'/','D_',FileAVI(1:end-4)],'D')
            if copyfile(FileAVI,FileAVI(1:end-4)) == 1
               delete(FileAVI); 
            end
        end
    catch
        %Code reaches this point if there is any errors during the video
        %processing and moves the error videos to a seperate folder
        if mkdir('ErrorVideos') == 1
            dest = 'ErrorVideos';
            if copyfile(FileAVI,dest) == 1
               delete(FileAVI); 
            end                        
        end
    end
end
toc