function invector = fourpamunmapA(outvector,mode)
    l = length(outvector);

    invector = '';
    if mode
        cutoff = [1.6 0.9 0];
    else
        cutoff = [4 2.5 1.1];
    end
    
    for k = 0:2:l*2-1

        if outvector(k/2 + 1) >= cutoff(1)
            invector(k+1) = '1';
            invector(k+2) = '1';
        elseif outvector(k/2 + 1) >= cutoff(2)
            invector(k+1) = '1';
            invector(k+2) = '0';
        elseif outvector(k/2 + 1) >= cutoff(3)
            invector(k+1) = '0';
            invector(k+2) = '1';
        else
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