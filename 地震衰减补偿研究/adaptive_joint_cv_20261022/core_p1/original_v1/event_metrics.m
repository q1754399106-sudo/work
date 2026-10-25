function T = event_metrics(D,names,c)
% Signed template amplitudes and local peak frequency; truth only for scoring.
% Windows are isolated here; not reliable for strongly overlapping thin beds.
rows=[];
for ix=1:numel(c.offset)
    for j=1:numel(c.event_t)
        keep=abs(c.t-c.event_t(j))<=0.07;
        u=c.t(keep)-c.event_t(j); z=pi*c.fpeak*u;
        w=(1-2*z.^2).*exp(-z.^2);
        truth=D(keep,ix,1); amp0=w'*truth/(w'*w);
        for k=[1 2 4 5 6 7]
            y=D(keep,ix,k); amp=w'*y/(w'*w);
            taper=0.5-0.5*cos(2*pi*(0:numel(y)-1)'/(numel(y)-1));
            spec=abs(fft(y.*taper,4096)); f=(0:2048)'/(4096*c.dt);
            band=f>=5 & f<=150; ff=f(band); ss=spec(1:2049); ss=ss(band);
            [~,im]=max(ss);
            rows=[rows;c.offset(ix),j,k,amp,amp0,abs(amp-amp0)/max(abs(amp0),eps),ff(im)]; %#ok<AGROW>
        end
    end
end
T=array2table(rows,'VariableNames',{'Offset_m','Event','MethodIndex','SignedTemplateAmplitude','TrueTemplateAmplitude','RelativeAmplitudeError','PeakFrequency_Hz'});
T.Method=string(names(T.MethodIndex))';
end
