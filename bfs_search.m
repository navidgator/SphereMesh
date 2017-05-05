clc, clear, close all;

quiver2 = @(P,V,varargin) quiver( P(1,:), P(2,:), V(1,:), V(2,:), varargin{:} );
plot2 = @(P,varargin) plot(P(1,:), P(2,:), varargin{:} );

I = imread('body1.jpg');
% I = imread('body2.jpeg');

I = flipud(I);

I = I(:,:,1) / 255; % Coverting image pixels to black and white
%imagesc(I);

S = double(I);       % Using for any function that needs "double image"    
S2 =  1-S;      % fliping the image color (black->white  , white -> black)

mex  CXXFLAGS='$CXXFLAGS -std=c++11' dtform.cpp;


[DIST, CORRS] = dtform(S);  % dtform for the Original image 

slope_x = conv2(DIST,[-1,0,1],'same'); % x-derivative of distance transform
slope_y = conv2(DIST,[-1;0;1],'same'); % y-derivative of distance transform
 
skeleton_x = -conv2(slope_x,[-1,0,1],'same'); % 2nd x-derivative of distance transform
skeleton_x(skeleton_x<0) = 0;

skeleton_y = -conv2(slope_y,[-1;0;1],'same'); % 2nd y-derivative of distance transform
skeleton_y(skeleton_y<0) = 0;

skeleton = DIST .* (skeleton_y.^2 + skeleton_x.^2); % skeleton of shape

skeleton(skeleton < 50) = 0;
skeleton(skeleton > 0) = 1;

[SKELDIST, SKELCORRS] = dtform(skeleton); % distance transform of skeleton

figure; axis equal; hold all; imagesc(SKELDIST);

figure; axis equal; hold all; imagesc(skeleton + S);

start_point = round(rand(1,2) .* size(I)); % random point to start search
[start_point(1),start_point(2)] = corrs_project(start_point,SKELCORRS);

skeleton_points = start_point; % centers of spheres in the sphere mesh
edges = []; % indices of connected spheres. each row defines a pill
frontier = [1;1]; % indices of spheres that have not been searched from yet

n = 1;
while size(frontier) ~= 0 % search until all possibile new points eliminated
    
    current_idx = frontier(1,1);
    parent_idx = frontier(2,1);
    
    edges = cat(1,edges,[parent_idx,current_idx]);
    
    frontier = frontier(:,2:end);
    
    current = skeleton_points(current_idx,:);
    
    % number of points on circle to start search
    search_n = 32;
    
    f = 1.2; % multiple of sphere radius to use as step size
    
    new_points = [];
    
    radius = DIST(current(1),current(2));
    if radius == 0 % do not search from zero-radius spheres
        continue;
    end
    
    theta = linspace(0,2*pi,search_n);
    for i = 1:search_n 
        
        new_point = round(current + f * radius * [cos(theta(i)),sin(theta(i))]);
        new_point(1) = max([min([new_point(1),size(I,1)]),1]);
        new_point(2) = max([min([new_point(2),size(I,2)]),1]);
            
        delta = inf;
        
        its = 0;
        % search for local maxima of distance transform
        while delta > 0.1
            
            offset = new_point - current;
            reprojected = current + f * radius * offset / norm(offset);
            
            reprojected = round(reprojected);
            reprojected(1) = max([min([reprojected(1),size(I,1)]),1]);
            reprojected(2) = max([min([reprojected(2),size(I,2)]),1]);
            
            [rx,ry] = corrs_project(reprojected,SKELCORRS);
            reprojected = [rx,ry];
            
            delta = norm(new_point - reprojected);
            new_point = reprojected;
            its = its + 1;
            
            if its > 50
                break;
            end
            
        end
        if delta < 1.5 && norm(new_point - current) < radius * f * 1.3
            new_points = cat(1,new_points,new_point);
        end
    end
    
    culled_points = [];
    k = 1;
    for i = 1:size(new_points,1)
        cull = false;
        for j = 1:size(culled_points,1)
            
            delta_theta = norm(new_points(i,:)-culled_points(j,:))/radius;
            
            if delta_theta < 0.5 % remove search points that converged to same angle
                cull = true;
                break;
            end
            
        end
        if ~cull
            for j = 1:size(skeleton_points,1)
                if j ~= current_idx
                    radius_j = DIST(skeleton_points(j,1),skeleton_points(j,2));
                    
                    % remove points that are too close to existing spheres
                    if norm(new_points(i,:) - skeleton_points(j,:)) < 0.5 * radius_j * f
                        cull = true;
                        break;
                    end
                end
            end
        end
        
        if ~cull
            culled_points(k,:) = new_points(i,:);
            k = k + 1;
        end
    end
    
    new_points = culled_points;
    
    skeleton_points(n+1:n+size(new_points,1),:) = new_points;
    
    % add new points to frontier list to be searched
    frontier = cat(2,frontier,[n+1:n+size(new_points,1);current_idx*ones(1,size(new_points,1))]);
    
    n = n + size(new_points,1);
    
end

figure; axis equal; hold all; imagesc(I);

scatter(skeleton_points(:,2),skeleton_points(:,1));

for i = 1:size(skeleton_points,1)
    
    draw_circle(skeleton_points(i,:),DIST(skeleton_points(i,1),skeleton_points(i,2)));
    
end

start_points = skeleton_points(edges(:,1),:);
end_points = skeleton_points(edges(:,2),:);

for i = 1:size(edges,1)
    
    c1 = start_points(i,:);
    c2 = end_points(i,:);
    
    plot([c1(2),c2(2)],[c1(1),c2(1)],'y');
    
    r1 = DIST(c1(1),c1(2));
    r2 = DIST(c2(1),c2(2));
    
    [R1,R2,L1,L2] = Bitangent (c1,r1,c2,r2);
    plot(R1(2),R1(1),'.R','MarkerSize',5);
    plot(R2(2),R2(1),'.R','MarkerSize',5);
    plot(L1(2),L1(1),'.R','MarkerSize',5);
    plot(L2(2),L2(1),'.R','MarkerSize',5);
    
    plot ([R1(2),R2(2)],[R1(1),R2(1)],'k');
    plot ([L1(2),L2(2)],[L1(1),L2(1)],'k');
    
end

function [xp,yp] = corrs_project(p,CORRS)

    [xp,yp] = ind2sub(size(CORRS),CORRS(p(1),p(2)));

end

function h = draw_circle(center, radius)
    
    th = 0:pi/100:2*pi;
    xunit = radius * cos(th) + center(2);
    yunit = radius * sin(th) + center(1);
    h = plot(xunit, yunit,'w');

end
