function plot_particles_2d(fig,L,r,x,y)
    % number of particles
    N       = length(r);

    % determine colors
    c       = zeros(N,3);
    c0      = [0 0.2 0.95];
    rmin    = min(r);
    rs      = r./rmin;        
    for n = 1:N
        c(n,:) = rs(n).*c0;
    end
    cmax    = max(max(c));
    c       = c./cmax; 

    figure(fig), clf, hold on, box on;
    axis('equal');
    axis([0 L(1) 0 L(2)]);

    x = mod(x,L(1));
    y = mod(y,L(2));
    for n = 1:N
        
        rectangle('Position',[x(n)-r(n), y(n)-r(n), 2*r(n), 2*r(n)],'Curvature',[1 1],'edgecolor',c(n,:),'facecolor',c(n,:));
        
        if (x(n)+r(n))>L(1)
            rectangle('Position',[x(n)-r(n)-L(1), y(n)-r(n), 2*r(n), 2*r(n)],'Curvature',[1 1],'edgecolor',c(n,:),'facecolor',c(n,:));
        end
        
        if (x(n)-r(n))<0
            rectangle('Position',[x(n)-r(n)+L(1), y(n)-r(n), 2*r(n), 2*r(n)],'Curvature',[1 1],'edgecolor',c(n,:),'facecolor',c(n,:));
        end
        
        if (y(n)+r(n))>L(2)
            rectangle('Position',[x(n)-r(n), y(n)-r(n)-L(1), 2*r(n), 2*r(n)],'Curvature',[1 1],'edgecolor',c(n,:),'facecolor',c(n,:));
        end
        
        if (y(n)-r(n))<0
            rectangle('Position',[x(n)-r(n), y(n)-r(n)+L(1), 2*r(n), 2*r(n)],'Curvature',[1 1],'edgecolor',c(n,:),'facecolor',c(n,:));
        end        
        
        if ((x(n)+r(n))>L(1) || (x(n)-r(n))<0) && ((y(n)+r(n))>L(2) || (y(n)-r(n))<0)
            x1 = mod((x(n) + L(1)),L(1));
            y1 = mod((y(n) + L(2)),L(2));
            rectangle('Position',[x1-r(n), y1-r(n), 2*r(n), 2*r(n)],'Curvature',[1 1],'edgecolor',c(n,:),'facecolor',c(n,:));
        end
        
    end
    drawnow;
    %hold off;
end












