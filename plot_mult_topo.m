function plot_mult_topo(mat,chans,desc,lims,time_vec,fname,fps)
% mat:                  Rows: Timeseries
%                       Cols: Channels
%                       If more than two dimensions are provided, topoplots
%                       will be plottet into subplots
% chans:                Channel Locations 
% desc:                 Description/ Title of the plot
% lims:                 Limits of the colorbar, either als two column
%                       vector or matrix if different limits should be
%                       applied for the different plots
%time_window and noverlap: settings of the calculations
% fname                 filename of the produced video 'Test.avi'
% fps:                  Frames per second, e.g. 10

                        
%% Define number of subplots
num_plots = size(mat,3);
rows = ceil(num_plots/2);
cols = ceil(num_plots/rows);
%% Check input arguments
if size(lims,1) == 1
    lims = repmat(lims,[num_plots,1]);
end
%% Define number of frames
framesNo = size(mat,1);
% allocate memory for future figures
F = struct('cdata', cell(1,framesNo), 'colormap', cell(1,framesNo));
pb = CmdLineProgressBar(("Preparing "+framesNo+ " Frames... "));

%% Create Figure and axes
figure('WindowState','maximized','Visible','on')
for i = 1:framesNo
    pb.print(i,framesNo)
    for j = 1:num_plots

        if num_plots == 1
            t1 = [num2str(time_vec(i)), ' s'];
            t2 = desc;
            t = [t1;t2];
            title(t);
        else
          ax(j) = subplot(rows,cols,j);
          sgtitle([num2str(time_vec(i)), ' s'])
          title(desc{1,j});
          
        end
        topoplot(mat(i,:,j),chans);
        colorbar
        caxis([lims(j,1),lims(j,2)])
        F(i) = getframe(gcf);
        if num_plots == 1
        cla(gca)
        end
    end
end

writerObj = VideoWriter(fname);
writerObj.FrameRate = fps;
disp("Creating " + length(F) / writerObj.FrameRate + " s video...");
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i);
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
disp('Done.')
disp(strcat('Find it as ', pwd, filesep, fname))
end