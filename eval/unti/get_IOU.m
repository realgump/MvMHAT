function o = get_IOU(B1, B2)
% Compute the symmetric intersection over union overlap between a set of
% bounding boxes in a and a single bounding box in b.
%
% a  a matrix where each row specifies a bounding box
% b  a matrix where each row specifies a bounding box

% AUTORIGHTS
% -------------------------------------------------------
% Copyright (C) 2011-2012 Ross Girshick
% Copyright (C) 2008, 2009, 2010 Pedro Felzenszwalb, Ross Girshick
% 
% This file is part of the voc-releaseX code
% (http://people.cs.uchicago.edu/~rbg/latent/)
% and is available under the terms of an MIT-like license
% provided in COPYING. Please retain this notice and
% COPYING if you use this file (or a portion of it) in
% your project.
% -------------------------------------------------------
    a = B1;
    a(3) = B1(1) + B1(3);
    a(4) = B1(2) + B1(4);
    b = B2;
    b(3) = B2(1) + B2(3);
    b(4) = B2(2) + B2(4);

    x1 = max(a(1), b(1));
    y1 = max(a(2), b(2));
    x2 = min(a(3), b(3));
    y2 = min(a(4), b(4));

    w = x2-x1+1;
    h = y2-y1+1;
    
    if w < 0 || h < 0
        o = 0;
    else
    
        inter = w.*h;
        aarea = (a(3)-a(1)+1) .* (a(4)-a(2)+1);
        barea = (b(3)-b(1)+1) * (b(4)-b(2)+1);
        % intersection over union overlap
        o = inter / (aarea+barea-inter);

    end
    
end
