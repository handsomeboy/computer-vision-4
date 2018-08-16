function [F0] = findFGroundTruth(MintL,MextL,MintR,MextR) 
% Input - Left and right intrinsic and extrinsic calibration matrices 
%          Mint and Mext
% Output - F0 - The ground truth fundamental matrix

PIR = MintR*MextR;
PIL = MintL*MextL;

m = null(PIR); %% Should be PI1 

PILm = PIL*m;

PILm_cp = [ 0, -PILm(3), PILm(2);
           PILm(3), 0, -PILm(1);
           -PILm(2), PILm(1), 0];
       
F0 = PILm_cp*PIL*PIR'/(PIR*PIR');

end