function img_proj = PlotCellsColorIndex(SFP, index)
%PLOTCELLSCOLORINDEX Summary of this function goes here
%   Detailed explanation goes here

for cell_i = 1:length(index)
    
    current_SFP = SFP(:,:,index(cell_i));
    current_SFP = current_SFP./max(current_SFP(:)); % Normalize
    
    current_color = [1-cell_i/length(index) 0 cell_i/length(index)];
    
    temp_img(:,:,1) = current_SFP*current_color(1); % Red channel
    temp_img(:,:,2) = current_SFP*current_color(2); % Green channel
    temp_img(:,:,3) = current_SFP*current_color(3); % Blue channel
    
    if cell_i == 1
        img_proj = temp_img;
    else
        img_proj = imadd(img_proj,temp_img);
    end

    imagesc(img_proj)
    drawnow
    pause
    
end

