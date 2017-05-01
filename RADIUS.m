function r =  RADIUS (N,p,p_i)
p12 = p_i - p;
len = norm(p12);
cost = dot(p12/ len , -N);
if (cost==0)
%     plot(p(2),p(1),'.r','MarkerSize',20);
%     plot(p_i(2),p_i(1),'.g','MarkerSize',20);
%     quiver(p(2),p(1),N(2),N(1),50);
r = 0;
return;
end
r = abs (len / (cost * 2));

end