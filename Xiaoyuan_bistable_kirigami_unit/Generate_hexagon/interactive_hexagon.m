% Create an interactive figure for deploying hexagon unit depending on
% deploying angle phi
function interactive_hexagon(edgeLen,l1,l2,l3,l4,beta,t)
    % Create a figure
    fig = figure('Position', [100, 100, 800, 800]); % Set figure size: left, bottom, width, height

    xlim([-l1*3, l1*3]); 
    ylim([-l1*3, l1*3]); 

    % Initial delta value
    delta = 0;
    prev_alpha_1 = pi/3;
    prev_alpha_2 = 2*pi/3;

    % Set the range of phi value
    strain_min = 0;
    strain_max = 0.8;
    
    % Define colors
    % colour = {'white', 'black', 'black', 'black'};
    colour = {'white', [206,101,95]/255, [90,174,52]/255, [109,131,250]/255};

    % Create the hexagon unit
    [hexagon,alpha_1_optimal,alpha_2_optimal] = hexagon_unit(prev_alpha_1, prev_alpha_2, delta, beta, edgeLen, l1,l4,t);

    % Add slider in the bottom half
    slider_handle = uicontrol(fig, 'Style', 'slider', 'Min', strain_min*edgeLen, 'Max', strain_max*edgeLen, 'Value', delta, ...
        'Position', [150, 30, 500, 30], 'Callback', {@update_plot, fig});
    % Add text label for slider value
    text_label = uicontrol(fig, 'Style', 'text', 'String', ['Strain = ', num2str(delta/edgeLen)], ...
        'Position', [650, 25, 60, 40],'FontSize', 12);

    % Update plot function
    function update_plot(hObject, eventdata, fig_handle)
        % Get the current phi value from the slider
        delta = get(hObject, 'Value');

        % Create the updated hexagon unit
        hexagon = hexagon_unit(alpha_1_optimal, alpha_2_optimal,delta, beta, edgeLen, l1,l4,t);

        % Clear the current plot and redraw
        clf(fig);
        axis equal;
        hold on;
        xlim([-l1*3, l1*3]);
        ylim([-l1*3, l1*3]);

        % Plot the updated hexagon unit
        plot_hexagon(hexagon,colour);

        % Add slider in the bottom half
        slider_handle = uicontrol(fig, 'Style', 'slider', 'Min', strain_min*edgeLen, 'Max', strain_max*edgeLen, 'Value', delta, ...
            'Position', [150, 30, 500, 30], 'Callback', {@update_plot, fig});

        % Add text label for slider value
        text_label = uicontrol(fig, 'Style', 'text', 'String', ['Strain = ', num2str(delta/edgeLen)], ...
            'Position', [650, 25, 60, 40],'FontSize', 12);
    end

    % Plot initial hexagon
    plot_hexagon(hexagon_unit(alpha_1_optimal, alpha_2_optimal,delta, beta, edgeLen, l1,l4,t),colour);
end


