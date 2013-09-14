
clear all;

MAKE_SPECKLE_IMAGES = false;
TIME_AVG_SPECKLE = false;


NUM_FILES =  207;

% The speckle pattern at t=0 of the simulation 
% (i.e. the speckle pattern without any influence of ultrasound).
t0 = dlmread('speckle_t0.dat');

% Data is in matrix of NxN (CCD grid dimensions).  Data is reshaped to fit in container.
CCD_xdim = size(t0,1);
CCD_ydim = size(t0,2);
speckle_data   = zeros(NUM_FILES, CCD_xdim*CCD_ydim);

delta_contrast = [];
delta_contrast_steffen = [];
avg = zeros(	CCD_xdim, CCD_ydim);

t0_1Darray = reshape(t0, 1, CCD_xdim*CCD_ydim);
for i=2:NUM_FILES
	tn = dlmread(['speckle_t', num2str(i-1), '.dat']);
	speckle_data(i,:) = reshape(tn, 1, size(tn,1)*size(tn,2));
	tn_1Darray = speckle_data(i,:);
	%delta_contrast = [delta_contrast; (std(tn)/mean(tn)- std(t0)/mean(t0))];
	%delta_contrast = [delta_contrast; abs(std(tn)/mean(tn)- std(t0)/mean(t0))];
	%delta_contrast = [delta_contrast; (std(t0)/mean(t0) - std(tn)/mean(tn))];
	
	% Calculate delta speckle contrast
	% ------------------------------------------------------------------------------
	delta_contrast = [delta_contrast; std(t0_1Darray)/mean(t0_1Darray) - ...
					       std(tn_1Darray)/mean(tn_1Darray)];

	% Calculate using Steffen's approach
	% ------------------------------------------------------------------------------
	%delta_contrast_steffen = [delta_contrast_steffen; mean((tn_1Darray - t0_1Darray).^2)/8];
	%if (i > 1)
		delta_contrast_steffen = [delta_contrast_steffen; mean((tn_1Darray/mean(tn_1Darray) - speckle_data(i-1,:)/mean(speckle_data(i-1,:))).^2)/8];	
	%end
	
	
	avg = tn + avg;
end

avg = avg ./ NUM_FILES;

%figure; plot(delta_contrast ./ max(delta_contrast), '-r');
%figure; plot(delta_contrast_steffen ./ max(delta_contrast_steffen), '-b');
%figure; plot(delta_contrast_steffen, '-b');





% Average speckle pattern over time (i.e. time_steps_avg)
% -------------------------------------------------------
if (TIME_AVG_SPECKLE)
	time_steps_avg = 5;
	avg_time_speckle = zeros(NUM_FILES/time_steps_avg, size(speckle_data,2));


	rows = size(speckle_data,1);
	A = mean(reshape(speckle_data, [rows/time_steps_avg, time_steps_avg, CCD_xdim*CCD_ydim]), 2);  % Perform the avg.
	A = reshape(A, [rows/time_steps_avg, CCD_xdim*CCD_ydim]);


	delta_contrast_avg_frames = zeros(1, size(A,1));
	size(delta_contrast_avg_frames)
	for i=1:size(A,1)
		delta_contrast_avg_frames(i) = std(reshape(A(1,:), CCD_xdim, CCD_ydim))/mean(reshape(A(1,:), CCD_xdim, CCD_ydim)) - ...
				      				   std(reshape(A(i,:), CCD_xdim, CCD_ydim))/mean(reshape(A(i,:), CCD_xdim, CCD_ydim)); 
	end
	figure; plot(delta_contrast_avg_frames);
end


% Form the speckle images and save them to a directory
% ----------------------------------------------------
if (MAKE_SPECKLE_IMAGES)
	set(0, 'defaultfigurevisible', 'off');  % disable displaying the plot.  Speeds things up considerably.
	figure;
	for i=0:NUM_FILES
		tn = dlmread(['speckle_t', num2str(i), '.dat']);
		imagesc(tn);
		filename=sprintf('images/%05d.png',i);
    	print(filename);
	end
end

% Two ways to make movie:
% 1) Mencoder
% mencoder mf://*.png -mf w=640:h=480:fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o output.avi
%
% 2) ffmpeg
% ffmpeg -i "%d5.png" -y output.mpeg
