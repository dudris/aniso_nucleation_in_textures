function [cmap , simtype_short] = plotMC2D_get_colormap_and_annotation(resentry)

    edge_intensity = 0.4;
    center_intensity = 0.9;
    single_color_intensity = edge_intensity*sqrt(3);

    if ~resentry.nucleation_on
        c_edge = [0, 0, single_color_intensity]; 
        simtype_short = 'NN';
    elseif resentry.nucleation_on && strcmp(resentry.nucl_mode,'uniform')
        c_edge = [single_color_intensity,0 , 0]; 
        simtype_short = 'IN';
    elseif resentry.nucleation_on && strcmp(resentry.nucl_mode,'aniso')
        c_edge = [0, single_color_intensity, 0]; 
        simtype_short = 'AN';
    end
    
%     disp(['edge BW intensity: ' num2str(sum(c_edge)/sqrt(3))]) % to display BW intensity
    c_center = center_intensity*[1, 1, 1];
    numcolors = 50;
    cmap = zeros(numcolors,3);
    % __ linear interpolation between c_edge an c_center and back
    for colmn = 1:3
        cmap(:,colmn) = [ linspace(c_edge(colmn),c_center(colmn),floor(numcolors/2))' ; linspace(c_center(colmn),c_edge(colmn),ceil(numcolors/2))' ];
    end
    
end % func