function plot_single_DSA(data,freq,xax)
    hold on
    pcolor(xax,freq,10*log10(data))
    shading flat
    shading interp
    colormap turbo
    colorbar
    caxis([-20 20])
    ylim([0 47])
    xlim([xax(1),xax(end)])
    ylabel('Frequency [Hz]')
    xlabel('Time [s]')
end