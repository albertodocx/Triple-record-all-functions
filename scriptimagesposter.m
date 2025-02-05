

figure(1)

subplot(3, 1, 1)
    imagesc(time_spec, imgfreqs, pow2db(thetach_profile(imrg, :)))
    title('SLM - Power Spectrum (dB)')
    ylabel('Fq (Hz)')
    xlabel('Time (s)')
%     colorbar
    colormap('jet')
    set(gca, 'YDir', 'normal')
    clim([0 100])

subplot(3, 1, 2)
    plot(time_spec, thetameansc, 'Color' ,'#dad0d8')
    hold on
    plot(time_spec, thetafilt, 'blue')
    title('SLM. Normalized Theta Power')
    xlim([0 time_spec(end)])
    xlabel('Time (s)')
    hold off

subplot(3, 1, 3)
    plot(time_spec, ratiosc, 'Color' ,'#dad0d8')
    hold on
    plot(time_spec, ratiofilt, 'blue')
    title('Normalized Ratio Theta/Delta Power')
    xlim([0 time_spec(end)])
    xlabel('Time (s)')
    yline(LIAthresh, 'red')
    yline(runthresh, 'green')
    hold off
