% function test_boundary()
% clc
% close all
% p		= [100, 200]';
% v		= [0.3, 1]';
% Imax	= [ 1200, 600];
% 
% [b] = farthest(p, v, Imax)
% 
% 
% end
% 
% function [b] = farthest(p, v, Imax)
% if norm(v)==0
% 	error('Error (1): The normal vector is zero!')
% 	return
% end
% 
% p=p(:);		% to make sure it is column vector
% v=v(:);		% to make sure it is column vector
% 
% Imin = [1;1];
% 
% corners = [ Imin(1),Imin(2)
% 	Imax(1),Imin(2)
% 	Imax(1),Imax(2)
% 	Imin(1),Imax(2)]';
% 
% end



function [b] = farthest(p, v, Imax)
if norm(v)==0
	error('Error (1): The normal vector is zero!')
	return
end

p=p(:);		% to make sure it is column vector
v=v(:);		% to make sure it is column vector

Imin = [1;1];

p_mat	= repmat(p,[1,4]);

corners = [ Imin(1),Imin(2)
	Imax(1),Imin(2)
	Imax(1),Imax(2)
	Imin(1),Imax(2)]';

d_mat	= corners - p_mat;

th	= atan2(v(1),v(2))*180/pi;
x	= atan2(d_mat(1,:),d_mat(2,:))*180/pi;
x	= x - th;
x	= x - sign(x).*floor(abs(x)./180).*360;
x	= [x 0];
[x,ind]	= sort(x);
[tmp,ind2] = ismember(5,ind);
if tmp==0
	error('Error (2)')
end



ind = [ind(ind2-1),ind(ind2+1)];
c_sel	= corners(:,ind);

d_sel		= c_sel(:,2) - c_sel(:,1);
[tmp,ind_xy] = ismember(0,d_sel);
if tmp==0
	error('Error (3)')
end

L	= (c_sel(ind_xy,2) - p(ind_xy))./v(ind_xy);

ind_yx	= 3 - ind_xy; 
b(ind_xy)	= c_sel(ind_xy,1);
b(ind_yx)	= p(ind_yx) + L*v(ind_yx);
b = b';

end
