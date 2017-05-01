clc, clear, close all;

quiver2 = @(P,V,varargin) quiver( P(1,:), P(2,:), V(1,:), V(2,:), varargin{:} );
plot2 = @(P,varargin) plot(P(1,:), P(2,:), varargin{:} );



I = imread('body1.jpg');
% I = imread('body2.jpeg');
% I = imread('fire.png'); 
% I = imread ('head.jpg');
% I = imread('girl2.jpg');

I = I(:,:,1) / 255; % Coverting image pixels to black and white
imagesc(I);

S = double(I);      % Using for any function that needs "double image"
S2 =  1-S;      % fliping the image color (black->white  , white -> black)

mex  CXXFLAGS='$CXXFLAGS -std=c++11' dtform.cpp;


[DIST, CORRS] = dtform(S);  % dtform for the Original image
figure ; imagesc(DIST);
[DIST2, CORRS2] = dtform(S2); % dtform for the color-fliped image
% figure; imagesc(DIST2);
figure; axis equal; hold all; 
[P,~] = contour(I, [0.5,0.5]);  % boundry of the image

%%--- Compute normals%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_shift = 10;
E = circshift(P,-n_shift,2)-P;
N1 = [0 -1;1 0]*E; % rotate edge by 90 degree
for i=1:size(N1,2), N1(:,i) = N1(:,i)./norm(N1(:,i)); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quiver2(P,N);

% r_ind = round (1 + size(P,2) * rand); % COMPUTING A RANDOM INDEX

% p_index = r_ind;        % THE INDEX OF THE FIRST POINT (Random)
% for p_index = 1 :size(P,2)

iter = round(size(P,2)/500);
% iter = 1;
figure; axis equal; hold all;
for p_index = 1 :iter:size(P,2)
%     p_index= 301;
    p = P(:,p_index);
    ttt= p;
    ttt(1) = p(2);
    ttt(2) = p(1);
    p = ttt;
    if (p(2)> size(I,2) || p(1)> size(I,1))
        continue;
    end
   
    N = N1(:,p_index);
    ttt= N;
    ttt(1) = N(2);
    ttt(2) = N(1);
    N = ttt;
    
    [center,rho,PI] = Circle(I, p' , N' ,CORRS, CORRS2,2);

    pos = [center(2)-rho center(1)-rho  2*rho 2*rho];
    r = rectangle('Position',pos,'Curvature',[1 1],'FaceColor',[0 .5 .5],'EdgeColor','none');
    
    
%     plot(center(2),center(1),'.g','MarkerSize',5);
%     th = 0:pi/500:2*pi;
%     xunit = rho * cos(th) + center(2);
%     yunit = rho * sin(th) + center(1);
% %     area(xunit,yunit);
%     h = plot(xunit, yunit);
%     plot(center(2),center(1),'.g','MarkerSize',10);
%     plot(PI(2),PI(1),'.black','MarkerSize',10);
%     hold on;
%     plot(p(2),p(1),'.r','MarkerSize',10);
    pause(0.1);
end




return ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



