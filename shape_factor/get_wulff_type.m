%% numerically determined values of Omg_cross
function wulff_type = get_wulff_type(n,Omg)
    
    types = {'weak','strong','strong earcross','general'};
    
    assert(Omg>=0)
    
    if Omg < 1
        wulff_type = types{1}; % weak
    else
        if n==2
            wulff_type = types{2}; % corners but no ear crossing
        elseif n==3
            wulff_type = types{2}; % corners but no ear crossing
            
        elseif n==4
            Omg_cross = 9; % soa when ears cross
            if Omg<Omg_cross % Omg smaller than Omg_Cross
                wulff_type = types{2}; % corners but no ear crossing
            else 
                wulff_type = types{3}; % corners AND ear crossing
            end % if Omg cross
        
        elseif n==6
            Omg_cross = 5.7477; % soa when ears cross
            if Omg<Omg_cross % Omg smaller than Omg_Cross
                wulff_type = types{2}; % corners but no ear crossing
            else 
                wulff_type = types{3}; % corners AND ear crossing
            end % if Omg cross
            
        end% if n
        
%         if datetime > datetime('03-Feb-2023 12:00:00')
%             wulff_type = types(4);
%         end
    end% if Omg weak
   
end