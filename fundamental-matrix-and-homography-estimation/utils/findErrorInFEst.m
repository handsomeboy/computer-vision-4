function [Ferror] = findErrorInFEst(p,region,F,F0)

count = 0;
nPts = size(p,2);

for i = 1:nPts
    if p(1,i) >= region(1) && p(1,i) <= region(3)
        if p(2,i) >= region(2) && p(2,i) <= region(4)
            
            count = count + 1;
            pCrop(:,count) = p(:,i);
            Epipolar0(count,:) = F0*pCrop(:,count);
            Epipolar(count,:) = F*pCrop(:,count);    
            
            temp = cropLineInBox(Epipolar0(count,1:2)',Epipolar0(count,3),region);

            if isnan(temp(1)) == 1
                Epipolar0(count,:) = [];
                Epipolar(count,:) = [];
                count = count - 1;
            else
                Lep(count,:) = temp(1,:); 
                Rep(count,:) = temp(2,:);
            end
        end
    end
end

% https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_an_equation
Ferror = zeros(count,1);
for j = 1:count
    a = Epipolar(j,1);
    b = Epipolar(j,2);
    c = Epipolar(j,3);
    d = sqrt(a^2 + b^2);
    x0L = Lep(j,1);
    y0L = Lep(j,2);
    x0R = Rep(j,1);
    y0R = Rep(j,2);
    
    perpDistL = abs(a*x0L + b*y0L + c) / d;
    perpDistR = abs(a*x0R + b*y0R + c) / d;
    
    Ferror(j) = max(perpDistL,perpDistR);
end

if count == 0
    Ferror = 0;
end

Ferror = max(Ferror);

end

    

            
            