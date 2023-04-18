%
%
% EMT 2017-03-20
% based on makePrettyAxis() but excludes all ticks that aren't limits or zero


function makeSimpleAxis(varargin)
    axh = gca;
    xOnly = false;
    yOnly = false;
    lineThickness = 1;
    fontSize = 20;
    assignargs(varargin);

    if isappdata(axh, 'drawAxisOrigLims')
        return;
    end

    drawX = ~yOnly;
    drawY = ~xOnly;

    if drawX
        xTick = get(axh, 'XTick');
        [xTick, xTickMask] = stripNonzeroMidpoints(xTick);
        xTickLabel = cellFromCharArray(get(axh, 'XTickLabel'));
        xTickLabel = xTickLabel(xTickMask);
        
        xLabel = get(get(axh, 'XLabel'), 'String');
    end

    if drawY
        yTick = get(axh, 'YTick');
        [yTick, yTickMask] = stripNonzeroMidpoints(yTick);
        yTickLabel = cellFromCharArray(get(axh, 'YTickLabel'));
        yTickLabel = yTickLabel(yTickMask);
        yLabel = get(get(axh, 'YLabel'), 'String');
    end
    
    axis(axh, 'off');
    box(axh, 'off');

    if drawX
        drawAxis(xTick, 'axisOrientation', 'h', 'tickLabels', xTickLabel, 'axisLabel', xLabel, ...
            'lineThickness', lineThickness,'fontSize', fontSize);
    end
    if drawY
        drawAxis(yTick, 'axisOrientation', 'v', 'tickLabels', yTickLabel, 'axisLabel', yLabel, ...
            'lineThickness', lineThickness,'fontSize', fontSize);
    end
end

function c = cellFromCharArray(ca)

    c = cell(size(ca,1), 1);
    for i = 1:size(ca,1)
        c{i} = strtrim(ca(i,:));
    end

end


function [ticksOut, tickMask] = stripNonzeroMidpoints(ticks)

tickMask = false(size(ticks));

tickMask(1) = true;
tickMask(end) = true;
tickMask(ticks==0) = true;

ticksOut = ticks(tickMask);

end
