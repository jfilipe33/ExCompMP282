%%% Arena and goal drawing
function [] = stencil()
    wx = [0.5 0.5 2 2 3.5 3.5];
    wy = [1.5 2.5 2.5 1.5 1.5 2.5];
    arenax = [0 0 4 4 0];
    arenay = [0 3 3 0 0];
    obs1x = [1 1 1.5 1.5 1];
    obs1y = [1 2 2 1 1];
    obs2x = [2.5 2.5 3 3 2.5];
    obs2y = [2 3 3 2 2];
    plot(arenax, arenay, 'b-', 'LineWidth', 3);
    hold on;
    plot(obs1x, obs1y, 'b-', 'LineWidth', 2);
    plot(obs2x, obs2y, 'b-', 'LineWidth', 2);
    plot(wx, wy, 'k--', 'LineWidth', 1.5);
    bound1 = viscircles([wx(2),wy(2)],0.1,'LineStyle',':');
    bound2 = viscircles([wx(3),wy(3)],0.1,'LineStyle',':');
    bound3 = viscircles([wx(4),wy(4)],0.1,'LineStyle',':');
    bound4 = viscircles([wx(5),wy(5)],0.1,'LineStyle',':');
    bound5 = viscircles([wx(6),wy(6)],0.1,'LineStyle',':');
    
