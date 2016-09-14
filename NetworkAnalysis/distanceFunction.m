function dist = distanceFunction(pos1, pos2, dist_type)
    if dist_type == 1
        dist = sqrt((pos2(1) - pos1(1))^2 + (pos2(2) - pos1(2))^2);
    elseif dist_type == 2
        dist = abs(pos2(1) - pos1(1)) + abs(pos2(2) - pos1(2));
    end
end