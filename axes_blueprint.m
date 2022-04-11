function axes_blueprint(patient, channels, tmin, tmax, isVisible)
    arguments
        patient
        channels
        tmin
        tmax
        isVisible {mustBeMember(isVisible, ["on", "off"])} = "on"
    end

fs = patient.fs;

assert((tmin*fs < tmax*fs) & (tmax*fs <= size(patient.data, 2)), ...
    "tmin and tmax should not exceed the signal duration!" + newline +...
    "Your signal duration: " + size(patient.data, 2)/fs+" s");
assert((tmin >=0) & (tmax > 0), "tmin and tmax must be positive!")

% hardcoded number of channels
xx = [patient.chanlocs(1:64).X];
yy = [patient.chanlocs(1:64).Y];

minval_x = min(xx);
minval_y = min(yy);

xx_offset = (xx + abs(minval_x)) / 250 + .07;
yy_offset = (yy + abs(minval_y)) / 250 + .15;

f = figure('WindowState', 'maximized', 'Color', 'black', 'Visible', isVisible);

% hardcoded number of channels
for i = [xx_offset;yy_offset;[1:64]]
    axes('Position',[i(1) i(2) .05 .05], ...
        'ButtonDownFcn', @lineCallback);
    
    plot(patient.data(i(3), tmin*fs+1:tmax*fs+1), 'color', '#ffa500', 'ButtonDownFcn', @lineCallback);
    title(channels{i(3)}, 'Color','white');
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'XColor', 'white')
    set(gca,'YColor', 'white')
    set(gca, 'Color', 'none')
end
