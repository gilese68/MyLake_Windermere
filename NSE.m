%-----------------------AgriMetSoft.com------------------------------------
 %    This is a function for calculate Nash–Sutcliffe model 
 %              efficiency coefficient
%
function nse = NSE( arrObs,arrPred )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
soor=sum((arrPred-arrObs).^2);
makh=sum((arrObs-mean(arrObs)).^2);
nse=1-(soor/makh);
end

