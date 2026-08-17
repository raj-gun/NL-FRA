function [freq_rng_nonlinear_ord,freq_rng_nonlinear_ord_NOFRF] = nonlinear_freq_range(N,a,b,len_fft,Fs)
% N - maximum nonlinear order, a - minimum frequency component, b - maximum frequency component
% calculates the output frequency range of a nonlinear system with order N
% when input has a frequency range [a,b], where a < b

if a > b
    disp('a should be greater than b');
    return;
end

freq_rng_nonlinear_ord = cell(N,1);
freq_rng_nonlinear_ord_NOFRF = cell(N,1);

% if a == b
%     a = 0;
%     ab_chk = 1;
% else
%     ab_chk = 0;
% end

for n = 1:N
    
    i = floor((n*a)/(a+b)) + 1;
    cond_i = (n*b)/(a+b);
    if cond_i >= i && n ~= 1 
        freq_rng_nonlinear_ord{n,1} = cell(1,i+1);
        for k = 0:i
            if k < i
                freq_rng_nonlinear_ord{n,1}{1,k+1} = [n*a - k*(a+b) , n*b - k*(a+b)];
            elseif k == i
                freq_rng_nonlinear_ord{n,1}{1,k+1} = [0 , n*b - i*(a+b)];
            end
        end
        
    else
        freq_rng_nonlinear_ord{n,1} = cell(1,i);
        for k = 0:i-1
            freq_rng_nonlinear_ord{n,1}{1,k+1} = [n*a - k*(a+b) , n*b - k*(a+b)];
        end
        
    end
    
%     if n == 1 && ab_chk == 1
%         freq_rng_nonlinear_ord{n,1}{1,1}(1) = b;
%     end
    
    %calculate the corresponding fft indexe ranges
    if len_fft ~=0 && Fs ~=0
        ind = length(freq_rng_nonlinear_ord{n,1});
        temp = cell(ind,1);
        for ind2 = 1:ind
            temp{ind2,1} = round(freq_rng_nonlinear_ord{n,1}{1,(ind - ind2 + 1)}(1) * (len_fft/Fs))+1 :1: round(freq_rng_nonlinear_ord{n,1}{1,(ind - ind2 + 1)}(2) * (len_fft/Fs))+1;
        end
        freq_rng_nonlinear_ord_NOFRF{n,1} = unique(cat(2,temp{:,1}));
        
    end
    
end



end