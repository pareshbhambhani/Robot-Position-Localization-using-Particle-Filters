%Function to output robot's Localization (Poses) using LIDAR readings and map data 
function myPose = particleLocalization(ranges, scanAngles, map, param)

% Number of poses to calculate
N = size(ranges, 2);

myPose = zeros(3, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Map Parameters 

% the number of grids for 1 meter.
myResol = param.resol;
% the origin of the map in pixels
myOrigin = param.origin; 

% The initial pose is given
myPose(:,1) = param.init_pose;

% Decide the number of particles, M.
M = 200;

%Cov of Gaussian noise in motion
Q = diag([0.0015,0.0015,0.0005]);

% Create M number of particles
P = repmat(myPose(:,1), [1, M]);

t = param.t; % Time vector
pose = param.pose; % Actual pose

w = ones(1,M)/M; % weights for each particle
correlation = zeros(1,M);

%Check video recording bool
if (param.captureBool)
    fileName = ['record_', datestr(now,'yyyy-mm-dd_HH-MM-SS'), '.mp4'];
    videoObj = VideoWriter(fileName,'MPEG-4');
    videoObj.Quality=100;
    open(videoObj);
end


%%%%%%%%%%%%%%%%% Figure Generation%%%%%%%%%%%%%%%
figure('units','pixels','position',get(groot,'ScreenSize'),'Name','Position Localization using Particle Filter','NumberTitle','off');

% Subplots for actual pose vs predicted pose
subplot(3,3,3);
grid; hold on; 
xActual = plot(t(1), pose(1,1),'k', 'LineWidth', 2);
xCalc = plot(t(1), myPose(1,1),'r', 'LineWidth', 1);
ylabel('$x~$', 'FontSize', 26, 'Interpreter', 'latex');
title('Pose comparison in local co-ordinate frame', 'FontSize', 26, 'Interpreter', 'latex');

h1 = legend('Actual pose', 'Estimated pose' );
set(h1,'FontSize',18);

subplot(3,3,6); grid;
hold on; yActual = plot(t(1), pose(2,1),'k', 'LineWidth', 2);
yCalc = plot(t(1), myPose(2,1),'r', 'LineWidth', 2);
ylabel('$y~$', 'FontSize', 26, 'Interpreter', 'latex');

subplot(3,3,9); grid;
hold on; thetaActual = plot(t(1), pose(3,1),'k', 'LineWidth', 2);
thetaCalc = plot(t(1), myPose(3,1),'r', 'LineWidth', 2);
ylabel('$\theta~$', 'FontSize', 26, 'Interpreter', 'latex');
xlabel('$time~(s)$', 'FontSize', 20, 'Interpreter', 'latex'); 
set(findobj('type','axes'),'fontsize',18);

%Subplot for LIDAR readings on map
subplot(3,3,[2,5,8]);
imagesc(map);
colormap('gray');
hold on;
axis equal;
xlim([0 size(map,2)]);
ylim([0 size(map,1)]);
title('LIDAR readings');
lidar_global(:,1) =  (ranges(:,1).*cos(scanAngles + myPose(3,1)) + myPose(1,1))*myResol + myOrigin(1);
lidar_global(:,2) = (-ranges(:,1).*sin(scanAngles + myPose(3,1)) + myPose(2,1))*myResol + myOrigin(2);
lidarPlot = plot(lidar_global(:,1), lidar_global(:,2), 'g.'); 
posPlot = plot(myPose(1,1)*param.resol+param.origin(1), ...
    myPose(2,1)*param.resol+param.origin(2), 'r.-');

%Subplot for particles' visualizations
subplot(3,3,[1,4,7]);
imagesc(map); hold on;

colormap('gray');
axis equal;
title('Particles');
particlesPlot = scatter(P(1,:)*myResol+myOrigin(1),P(2,:)*myResol+myOrigin(2), 6, 'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5);



for j = 2:N 

    %Propagate the particles. Encoder measurements are not provided.
    %Instead using a distribution with zero mean and Q cov
     P = diag(myPose(:,j-1))*ones(3,M) +  mvnrnd([0;0;0],Q,M)';
     
     
     %For each particle in P,do
    for p = 1:M
        
        % Find occupied-measurement cells from LIDAR readings for current particle
        occCells = grids(ranges(:,j),scanAngles,P(:,p),myOrigin,myResol);
        del_occ =  occCells(:,1)<1 | occCells(:,2)<1 |  occCells(:,1) > size(map,2) |  occCells(:,2) > size(map,1);

        occCells(del_occ,:) = [];
        
        % Find the correlation index between occupied values from LIDAR and
        % occupied values from MAP
        occ_index = sub2ind(size(map),occCells(:,2),occCells(:,1)); %indicies of occ cells in single value index format
        occ_values = map(occ_index); %MAP values of occ cells
        correlation(1,p) = sum(occ_values(occ_values>=0.5)*10) + sum(occ_values(occ_values<=-0.2)*2); %correlation for each particle
    end

    % Update Particle weights based on correlation score
    correlation(correlation<0)= 0;
    w = correlation;
    w = w./sum(w); % Normalize weights
    
    % Choose best particle to update the pose
    [~,index]=max(w); %Find the index of max weight value
    myPose(:,j) = P(:,index);
    


    % Resample if the effective number of particles is smaller than a threshold
   Neff = (sum(w))^2 / sum(w.^2); %Nr will be 1 since w is normalized
    
    if Neff < 0.1*M
        %Using Inverse Transform Sampling
        %(http://www.cse.psu.edu/~rtc12/CSE598C/samplingSlides.pdf)
        edges = min([0 cumsum(w)],1); % Create CDF for weights and add 0 at beginning
        edges(end) = 1;                 % get the upper edge exact
        
        u1 = rand/M; % Random number from where to start counting evenly spaced 200 readings
        [~, idx] = histcounts(u1:1/M:1,edges); %Find the index of the edges in which u lies
        P = P(:,idx);                % extract new particles
        w = w(:,idx);               % Update weights
        w = w./sum(w);              % Normalize weights
    end



    % Update the figures
    particlesPlot.XData = P(1,:)*myResol+myOrigin(1);
    particlesPlot.YData = P(2,:)*myResol+myOrigin(2);
    
    lidarPlot.XData = (ranges(:,j).*cos(scanAngles + myPose(3,j)) + myPose(1,j))*myResol + myOrigin(1);
    lidarPlot.YData = (-ranges(:,j).*sin(scanAngles + myPose(3,j)) + myPose(2,j))*myResol + myOrigin(2);
    
    currentPoseX = myPose(1,j)*param.resol+param.origin(1);
    currentPoseY = myPose(2,j)*param.resol+param.origin(2);
    posPlot.XData = [posPlot.XData currentPoseX];
    posPlot.YData = [posPlot.YData currentPoseY];
    
    xActual.XData = t(1:j)'; xActual.YData = [xActual.YData pose(1,j)];
    yActual.XData = t(1:j)'; yActual.YData = [yActual.YData pose(2,j)];
    thetaActual.XData = t(1:j)'; thetaActual.YData = [thetaActual.YData pose(3,j)];
    xCalc.XData = t(1:j)'; xCalc.YData = [xCalc.YData myPose(1,j)];
    yCalc.XData = t(1:j)'; yCalc.YData = [yCalc.YData myPose(2,j)];
    thetaCalc.XData = t(1:j)'; thetaCalc.YData = [thetaCalc.YData myPose(3,j)];
    
    subplot(3,3,[1,4,7]);
    xlim([currentPoseX-40 currentPoseX+40])
    ylim([currentPoseY-40 currentPoseY+40])
    
    drawnow;

    %Capture Video
    if (param.captureBool)
        writeVideo(videoObj, getframe(gcf));
    end
    
    
end
if (param.captureBool)
    close(videoObj);

end

