function site_xy = nucleation(nucl_probability,randnr, spec,tres)
    
%     figure(98)
%     plot(nucl_probability,'ro')

    %__ with no nucleation barrier, uniform probability for all these oris     
    no_nucl_barr = nucl_probability>=tres;
    if any(no_nucl_barr) 
        oris_no_barr = find(no_nucl_barr);
        ind_nucltd = ceil(randnr*length(oris_no_barr));
        site_xy = oris_no_barr(ind_nucltd);
        return
    end

    % orientation from prob. distrib assumes that oris=res.bori=res.tori
    is3Dnucleation = ~isnan(nucl_probability);
    nucleated = nucl_probability>=randnr;

    if any(is3Dnucleation & nucleated)
%             possible_oris = oris(is3Dnucleation & nucleated );
        possible_oris = find(is3Dnucleation & nucleated ); % not orientation but indices
        narrowed_prob = nucl_probability(is3Dnucleation & nucleated );

        if length(possible_oris)==1
            site_xy = possible_oris;
        else
            indvec = 1:length(possible_oris);
            loc_ori = rand_tori(indvec,narrowed_prob);
            ind_sampled =  abs(indvec-loc_ori)==min(abs(indvec-loc_ori)); % should be unambiguous as ori is monotonic
%             ind_sampled =  abs(possible_oris-loc_ori)==min(abs(possible_oris-loc_ori)); % should be unambiguous as ori is monotonic
            site_xy = possible_oris(ind_sampled);
        end
    else %nothing happened
        site_xy = 0;
    end % 
    
end % function