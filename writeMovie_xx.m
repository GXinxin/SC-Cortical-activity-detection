function writeMovie_xx(M, filename, useFFmpeg)
%writeMovie - Make avi movie
%PURPOSE -- Make movie video from a colormapped matlab movie structure. Called by Iarr2avi.m and others and used in conjunction with output from timeColorMapProj.m
%USAGE -- 	writeMovie(M, 'filename.avi');
% M - A matlab specific 'movie' data structure
% filename - string, the 'filename.avi' that the data come from and the string from which the output filename will be formatted
%
% See also timeColorMapProj.m, Iarr2montage.m, myMovie2avi.m, Iarr2avi.m
%
%James B. Ackman 2014-12-31 10:46:39
% Xinxin Ge 06/22/16

if nargin < 3 || isempty(useFFmpeg), useFFmpeg = 1; end
if nargin < 2 || isempty(filename), filename = ['movie' datestr(now,'yyyymmdd-HHMMSS') '.avi']; end

disp(['Making ' filename '-----------'])


    if useFFmpeg
            rng('shuffle')
            tmpPath = ['wbDXtmp' num2str(round(rand(1)*1e09))];
            mkdir(tmpPath)

        szZ = numel(M);
        for fr = 1:szZ; %option:parfor
            tmpFilename = fullfile(tmpPath, sprintf('img%05d.jpg',fr));
            if isempty(M(fr).colormap)
                imwrite(M(fr).cdata,tmpFilename)
            else
                imwrite(M(fr).cdata,M(fr).colormap,tmpFilename)
            end
        end

        
        filePath = [pwd, '\', tmpPath];
        tic
        disp('ffmpeg running...')
        try
            %System cmd to ffmpeg:
            system(['ffmpeg -f image2 -i ' filePath filesep 'img%05d.jpg -vcodec mjpeg ' filename])
            %The call to ffmpeg can be modified to write something other than a motion jpeg avi video:
            %system('ffmpeg -f image2 -i img%05d.png a.mpg')
            rmdir(tmpPath,'s');
        catch
            rmdir(tmpPath,'s');
            error(errstr);
        end
        toc

    else
        tic
        disp('using video obj...')
        vidObj = VideoWriter(filename);
        open(vidObj);
        for i =1:numel(M)
            writeVideo(vidObj,M(i));
        end
        close(vidObj);
        toc
    end
