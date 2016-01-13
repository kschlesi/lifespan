function offset = offset(res)
% set ID offset for left hemisphere Lausanne regions
    switch res
        case 33,  offset = 41;
        case 60,  offset = 64;
        case 125, offset = 115;
    end
end
