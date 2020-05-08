function [rate] = bitter_rate(CNR,ENdB,ENdB_lim,fb,fb_lim,hard)

%If ENdB is acceptable, data rate is equal to the maximum
if ENdB >= ENdB_lim  
    rate = fb;
end

%If ENdB is too low, reduce the symbol rate until the ENdB is acceptable or
%until the minimum rate of the hardware is reached
if ENdB < ENdB_lim
    while fb > fb_lim && ENdB < ENdB_lim   
        fb = round(fb/2);
        EbNo = CNR * hard.BW / fb;
        ENdB = 10*log10(EbNo);
    end

    
    %NEED TO FIX THIS DEFINTION OF fb
    if fb < 1000 * log2(hard.M) * hard.R    %[ksps] Check if loop was broken because hardware limit was exceeded
        rate = 0;
    else            %Otherwise, used the reduced bit-rate
        rate = fb; 
    end
end

end
            