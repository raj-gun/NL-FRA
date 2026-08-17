function [bicx,bicx_diag_mean,Bspecx,Bspecx_diag_mean,waxis_diag] = BISx(u,y,Fs,nfft,wind,nsamp,overlap,c,displ)
[bicx,~] = bicoherx (u,y,y,nfft,wind,nsamp,overlap);
bicx_diag_mean = diag_mean(flipud(bicx));

[Bspecx,~] = bispecdx (u,y,y,nfft,wind,nsamp,overlap);
Bspecx_diag_mean = diag_mean(flipud(abs(Bspecx)));

[bic,~] = bicoher (y,nfft,wind,nsamp,overlap);
bic_diag_mean = diag_mean(flipud(bic));

[Bspec,waxis] = bispecd (y,nfft,wind,nsamp,overlap);
Bspec_diag_mean = diag_mean(flipud(abs(Bspec)));
waxis_diag = (-1:length(waxis)-2).*mean(diff(waxis)).*Fs;
if displ==1
    figure;
    subplot(2,2,1);surf(waxis.*Fs,waxis.*Fs,bicx);title('BICx');colorbar; shading interp;view ([0 0 90]);xlim([-20 20]);ylim([-20 20]);colormap('parula')
    subplot(2,2,2);surf(waxis.*Fs,waxis.*Fs,abs(Bspecx));title('BISx');colorbar;shading interp;view ([0 0 90]);xlim([-20 20]);ylim([-20 20]);colormap('parula')
    subplot(2,2,3);surf(waxis.*Fs,waxis.*Fs,bic);title('BIC');colorbar;shading interp;view ([0 0 90]);xlim([-20 20]);ylim([-20 20]);colormap('parula')
    subplot(2,2,4);surf(waxis.*Fs,waxis.*Fs,abs(Bspec));title('BIS');colorbar;shading interp;view ([0 0 90]);xlim([-20 20]);ylim([-20 20]);colormap('parula')
    
    figure;
    subplot(2,2,1);plot(waxis_diag,bicx_diag_mean,'Color',c);title('BICx');%axis([0, inf, -inf, inf]);
    subplot(2,2,2);plot(waxis_diag,Bspecx_diag_mean,'Color',c);%axis([0, N*f1, -inf, inf]);%set(gca, 'YScale', 'log');
    title('BISx');
    subplot(2,2,3);plot(waxis_diag,bic_diag_mean,'Color',c);title('BIC');%axis([0, inf, -inf, inf]);
    subplot(2,2,4);plot(waxis_diag,Bspec_diag_mean,'Color',c);%axis([0, N*f1, -inf, inf]);%set(gca, 'YScale', 'log');
    title('BIS');
end
end