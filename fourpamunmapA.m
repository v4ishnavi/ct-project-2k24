function invector = fourpamunmapA(outvector, mode)
    l = length(outvector);

    invector = '';
    if mode
        cutoff = [1.6 0.9 0];
    else
        cutoff = [4 2.5 1.1];
    end
    cutoff = iterative_threshold(outvector, 4);
    disp(cutoff)
    cutoff = sort(cutoff); 
    for k = 0:2:l*2-1
        val = zeros(1, 4);
        for i = 1:4
        val(i) = abs(outvector(k/2 + 1) - cutoff(i));
        end
        if min(val) == val(4)
            invector(k+1) = '1';
            invector(k+2) = '1';
        elseif min(val) == val(3)
            invector(k+1) = '1';
            invector(k+2) = '0';
        elseif min(val) == val(2)
            invector(k+1) = '0';
            invector(k+2) = '1';
        elseif min(val) == val(1)
            invector(k+1) = '0';
            invector(k+2) = '0';
        end
        % if 0 >= outvector(k/2 + 1) && outvector(k/2 + 1) > -2
        %     invector(k+1) = '0';
        %     invector(k+2) = '0';
        % elseif outvector(k/2 + 1) > 0 && outvector(k/2 +1) <= 2
        %     invector(k+1) = '0';
        %     invector(k+2) = '1';
        % elseif outvector(k/2 + 1) <= -2 
        %     invector(k+1) = '1';
        %     invector(k+2) = '0';
        % else
        %     invector(k+1) = '1';
        %     invector(k+2) = '1';
        % end
    end

end