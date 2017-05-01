function [p,N] = New_Point(DIST,CORRS,S,RS,RD,epsilon)


PX = RS;
distance = 0;
while(1)
    PX_Next = PX + RD' * distance;
    rounded = round(PX_Next);
    distance = DIST(rounded(1),rounded(2));
%     plot(rounded(2),rounded(1),'.g','MarkerSize',20);

    if ( distance <  0.5)
%         PX = PX - RD' * 2.0;
        roundedP = round(PX);
        [JJ,II] = ind2sub(size(S),CORRS(roundedP(1),roundedP(2)));
        JJ = JJ+1;
        p = [JJ,II];
        N =  [JJ,II] - PX ;
        N = N/norm(N);
        return;
    end
    
      PX = PX_Next; 
    
    
end

end