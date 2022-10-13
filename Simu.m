size = 3;
L = 300;
for x = 0:300/size:300
    line([x x], [0 300], [20 20], 'color', 'b');
end
for y = 0:300/size:300
    line([0 300], [y y], [20 20], 'color', 'b');
end
draw_X(1, 2, L, size);
draw_O(1, 3, L, size);
draw_X(2, 1, L, size);
draw_O(2, 2, L, size);
draw_X(2, 3, L, size);
draw_O(3, 1, L, size);
draw_X(3, 2, L, size);
draw_O(3, 3, L, size);
draw_O(1, 1, L, size);

function [t1, t2] = draw_link(X, Y, Z)
    l2=250; l3=200; l4=120; X0=-50; Y0=150; Z0=0;
    c2 = (X^2 + Y^2 - l2^2 - l3^2)/(2*l2*l3);
    s2 = sqrt(abs(1-c2^2));
    t2 = atan2(s2, c2);

    l1 = Z + l4;

    c1 = (l3*c2+l2)*Y - X*l3*s2;
    s1 = -(l3*c2+l2)*X - Y*l3*s2;
    t1 = atan2(s1, c1);

    X1 = X0; Y1 = Y0; Z1 = Z0 + l1;
    X2 = X1 - sin(t1)*l2; Y2 = Y1 + cos(t1)*l2; Z2 = Z1;
    X3 = X2-sin(t1+t2)*l3; Y3 = Y2 + cos(t1+t2)*l3 ; Z3 = Z1;
    X4 = X3; Y4 = Y3; Z4 = Z3-l4;
    
    Px = -l3*sin(t1+t2)-l2*sin(t1)-50;
    Py = l3*cos(t1+t2)+l2*cos(t1)+150;
    Pz = l1-l4;

    grid on
    xlabel('x')
    ylabel('y')
    zlabel('z')
    link1 = line([X0, X1], [Y0, Y1], [Z0, Z1], 'Color', "black", 'LineWidth', 2);
    hold on;
    link2 = line([X1, X2], [Y1, Y2], [Z1, Z2], 'Color', "blue", 'LineWidth', 2);
    hold on;
    link3 = line([X2, X3], [Y2, Y3], [Z2, Z3], 'Color', "yellow", 'LineWidth', 2);
    hold on;
    link4 = line([X3, X4], [Y3, Y4], [Z3, Z4], 'Color', "red", 'LineWidth', 2);
    hold on;
    plot3(Px, Py, Pz, ".", 'Color', "red");
    hold on;
    axis([-600 400 -400 600 0 300]) 
    pause(0.05);
    delete(link1);
    delete(link2);
    delete(link3);
    delete(link4);

end
function draw_O(i,j,L,size)
    r = (1/3) * (L/size);
    X_center =  50+(i-0.5)*(L/size);
    Y_center = -150+(j-0.5)*(L/size);
    X = X_center + r;
    Y = Y_center ;
    l2=250; l3=200;
    c2 = (X^2 + Y^2 - l2^2 - l3^2)/(2*l2*l3);
    s2 = sqrt(abs(1-c2^2));
    t2 = atan2(s2, c2);
    c1 = (l3*c2+l2)*Y - X*l3*s2;
    s1 = -(l3*c2+l2)*X - Y*l3*s2;
    t1 = atan2(s1, c1);
    di_ra(t1,t2);
    xuong(X_center + r*cos(0),Y_center + r*sin(0), 30)
    for t=0:0.1:2*pi
        X = X_center + r*cos(t);
        Y = Y_center + r*sin(t);
        [t1, t2] = draw_link(X, Y, 20);
    end
    len(X_center + r*cos(0),Y_center + r*sin(0), 20);
    di_ve(t1,t2);
    
    
end
function draw_X(i,j,L,size)
    r = (1/3) * (L/size);
    X_center =  50+(i-0.5)*(L/size);
    Y_center = -150+(j-0.5)*(L/size);
    X = X_center + r;
    Y = Y_center - r ;
    l2=250; l3=200;
    c2 = (X^2 + Y^2 - l2^2 - l3^2)/(2*l2*l3);
    s2 = sqrt(abs(1-c2^2));
    t2 = atan2(s2, c2);
    c1 = (l3*c2+l2)*Y - X*l3*s2;
    s1 = -(l3*c2+l2)*X - Y*l3*s2;
    t1 = atan2(s1, c1);
    di_ra(t1,t2);
    xuong(X_center + r,Y_center - r, 30);
    for t=-r:1:r
        [t1tr, t2tr] = draw_link(X_center-t, Y_center+t, 20);
    end
    len(X_center -r,Y_center + r, 20);
    X = X_center + r;
    Y = Y_center - r ;
    c2 = (X^2 + Y^2 - l2^2 - l3^2)/(2*l2*l3);
    s2 = sqrt(abs(1-c2^2));
    t2s = atan2(s2, c2);
    c1 = (l3*c2+l2)*Y - X*l3*s2;
    s1 = -(l3*c2+l2)*X - Y*l3*s2;
    t1s = atan2(s1, c1);
    move_X(t1tr, t1s, t2tr, t2s);
    xuong(X_center+r ,Y_center + r, 30);
    for t=-r:1:r
        [t1, t2] = draw_link(X_center-t, Y_center-t, 20);
    end
    len(X_center -r,Y_center -r, 20);
    di_ve(t1,t2);
end
function di_ve(t1_ht,t2_ht)
    l2=250; l3=200; l4=120; X0=-50; Y0=150; Z0=0;
    for t2=t2_ht:-0.05:0
        t1 = t1_ht;
        l1 = l4+30;

        X1 = X0; Y1 = Y0; Z1 = Z0 + l1;
        X2 = X0 - sin(t1)*l2; Y2 = Y0 + cos(t1)*l2; Z2 = Z1;
        X3 = X2-sin(t1+t2)*l3; Y3 = Y2 + cos(t1+t2)*l3 ; Z3 = Z1;
        X4 = X3; Y4 = Y3; Z4 = Z3-l4;

        grid on
        d1 = line([X0, X1], [Y0, Y1], [Z0, Z1], 'Color', "black", 'LineWidth',2);
        hold on;
        d2 = line([X1, X2], [Y1, Y2], [Z1, Z2], 'Color', "blue", 'LineWidth',2);
        hold on;
        d3 = line([X2, X3], [Y2, Y3], [Z2, Z3], 'Color', "yellow", 'LineWidth',2);
        hold on;
        d4 = line([X3, X4], [Y3, Y4], [Z3, Z4], 'Color', "red", 'LineWidth',2);
        hold on;
        axis([-600 400 -400 600 0 300]) 
        pause(0.05);
        delete(d1); delete(d2); delete(d3); delete(d4);
    end
    for t1=t1_ht:0.05:0
        t2 = 0;
        l1 = l4+30;

        X1 = X0; Y1 = Y0; Z1 = Z0 + l1;
        X2 = X0 - sin(t1)*l2; Y2 = Y0 + cos(t1)*l2; Z2 = Z1;
        X3 = X2-sin(t1+t2)*l3; Y3 = Y2 + cos(t1+t2)*l3 ; Z3 = Z1;
        X4 = X3; Y4 = Y3; Z4 = Z3-l4;

        grid on
        d1 = line([X0, X1], [Y0, Y1], [Z0, Z1], 'Color', "black", 'LineWidth',2);
        hold on;
        d2 = line([X1, X2], [Y1, Y2], [Z1, Z2], 'Color', "blue", 'LineWidth',2);
        hold on;
        d3 = line([X2, X3], [Y2, Y3], [Z2, Z3], 'Color', "yellow", 'LineWidth',2);
        hold on;
        d4 = line([X3, X4], [Y3, Y4], [Z3, Z4], 'Color', "red", 'LineWidth',2);
        hold on;
        axis([-600 400 -400 600 0 300]) 
        pause(0.05);
        delete(d1); delete(d2); delete(d3); delete(d4);
    end
end
function di_ra(t1_ht,t2_ht)
        l2=250; l3=200; l4=120; X0=-50; Y0=150; Z0=0;
        for t2=0:0.05:t2_ht
            t1 = 0;
            l1 = l4+30;

            X1 = X0; Y1 = Y0; Z1 = Z0 + l1;
            X2 = X0 - sin(t1)*l2; Y2 = Y0 + cos(t1)*l2; Z2 = Z1;
            X3 = X2-sin(t1+t2)*l3; Y3 = Y2 + cos(t1+t2)*l3 ; Z3 = Z1;
            X4 = X3; Y4 = Y3; Z4 = Z3-l4;

            grid on
            d1 = line([X0, X1], [Y0, Y1], [Z0, Z1], 'Color', "black", 'LineWidth',2);
            hold on;
            d2 = line([X1, X2], [Y1, Y2], [Z1, Z2], 'Color', "blue", 'LineWidth',2);
            hold on;
            d3 = line([X2, X3], [Y2, Y3], [Z2, Z3], 'Color', "yellow", 'LineWidth',2);
            hold on;
            d4 = line([X3, X4], [Y3, Y4], [Z3, Z4], 'Color', "red", 'LineWidth',2);
            hold on;
            axis([-600 400 -400 600 0 300]) 
            pause(0.05);
            delete(d1); delete(d2); delete(d3); delete(d4);
        end
        for t1=0:-0.05:t1_ht
            t2 = t2_ht;
            l1 = l4+30;

            X1 = X0; Y1 = Y0; Z1 = Z0 + l1;
            X2 = X0 - sin(t1)*l2; Y2 = Y0 + cos(t1)*l2; Z2 = Z1;
            X3 = X2-sin(t1+t2)*l3; Y3 = Y2 + cos(t1+t2)*l3 ; Z3 = Z1;
            X4 = X3; Y4 = Y3; Z4 = Z3-l4;

            grid on
            d1 = line([X0, X1], [Y0, Y1], [Z0, Z1], 'Color', "black", 'LineWidth',2);
            hold on;
            d2 = line([X1, X2], [Y1, Y2], [Z1, Z2], 'Color', "blue", 'LineWidth',2);
            hold on;
            d3 = line([X2, X3], [Y2, Y3], [Z2, Z3], 'Color', "yellow", 'LineWidth',2);
            hold on;
            d4 = line([X3, X4], [Y3, Y4], [Z3, Z4], 'Color', "red", 'LineWidth',2);
            hold on;
            axis([-600 400 -400 600 0 300]) 
            pause(0.05);
            delete(d1); delete(d2); delete(d3); delete(d4);
        end
end
function len(X, Y, Z)
    l2=250; l3=200; l4=120; X0=-50; Y0=150; Z0=0;
    for Zadd = 1:1:10
        
        c2 = (X^2 + Y^2 - l2^2 - l3^2)/(2*l2*l3);
        s2 = sqrt(abs(1-c2^2));
        t2 = atan2(s2, c2);

        l1 = Zadd + l4 + Z;

        c1 = (l3*c2+l2)*Y - X*l3*s2;
        s1 = -(l3*c2+l2)*X - Y*l3*s2;
        t1 = atan2(s1, c1);

        X1 = X0; Y1 = Y0; Z1 = Z0 + l1;
        X2 = X1 - sin(t1)*l2; Y2 = Y1 + cos(t1)*l2; Z2 = Z1;
        X3 = X2-sin(t1+t2)*l3; Y3 = Y2 + cos(t1+t2)*l3 ; Z3 = Z1;
        X4 = X3; Y4 = Y3; Z4 = Z3-l4;

        grid on
        xlabel('x')
        ylabel('y')
        zlabel('z')
        link1 = line([X0, X1], [Y0, Y1], [Z0, Z1], 'Color', "black", 'LineWidth', 2);
        hold on;
        link2 = line([X1, X2], [Y1, Y2], [Z1, Z2], 'Color', "blue", 'LineWidth', 2);
        hold on;
        link3 = line([X2, X3], [Y2, Y3], [Z2, Z3], 'Color', "yellow", 'LineWidth', 2);
        hold on;
        link4 = line([X3, X4], [Y3, Y4], [Z3, Z4], 'Color', "red", 'LineWidth', 2);
        hold on;
        axis([-600 400 -400 600 0 300]) 
        pause(0.05);
        delete(link1);
        delete(link2);
        delete(link3);
        delete(link4);
    end
end
function xuong(X, Y, Z)
    l2=250; l3=200; l4=120; X0=-50; Y0=150; Z0=0;
    for Zxuong = 1:1:10
        c2 = (X^2 + Y^2 - l2^2 - l3^2)/(2*l2*l3);
        s2 = sqrt(abs(1-c2^2));
        t2 = atan2(s2, c2);

        l1 = l4 + Z -Zxuong;

        c1 = (l3*c2+l2)*Y - X*l3*s2;
        s1 = -(l3*c2+l2)*X - Y*l3*s2;
        t1 = atan2(s1, c1);

        X1 = X0; Y1 = Y0; Z1 = Z0 + l1;
        X2 = X1 - sin(t1)*l2; Y2 = Y1 + cos(t1)*l2; Z2 = Z1;
        X3 = X2-sin(t1+t2)*l3; Y3 = Y2 + cos(t1+t2)*l3 ; Z3 = Z1;
        X4 = X3; Y4 = Y3; Z4 = Z3-l4;

        grid on
        xlabel('x')
        ylabel('y')
        zlabel('z')
        link1 = line([X0, X1], [Y0, Y1], [Z0, Z1], 'Color', "black", 'LineWidth', 2);
        hold on;
        link2 = line([X1, X2], [Y1, Y2], [Z1, Z2], 'Color', "blue", 'LineWidth', 2);
        hold on;
        link3 = line([X2, X3], [Y2, Y3], [Z2, Z3], 'Color', "yellow", 'LineWidth', 2);
        hold on;
        link4 = line([X3, X4], [Y3, Y4], [Z3, Z4], 'Color', "red", 'LineWidth', 2);
        hold on;
        axis([-600 400 -400 600 0 300]) 
        pause(0.05);
        delete(link1);
        delete(link2);
        delete(link3);
        delete(link4);
    end
end
function move_X(t1tr, t1s, t2tr, t2s)
l2=250; l3=200; l4=120; X0=-50; Y0=150; Z0=0;
    for t2=t2tr:-0.05:t2s
        t1 = t1tr;
        l1 = l4+30;

        X1 = X0; Y1 = Y0; Z1 = Z0 + l1;
        X2 = X0 - sin(t1)*l2; Y2 = Y0 + cos(t1)*l2; Z2 = Z1;
        X3 = X2-sin(t1+t2)*l3; Y3 = Y2 + cos(t1+t2)*l3 ; Z3 = Z1;
        X4 = X3; Y4 = Y3; Z4 = Z3-l4;

        grid on
        d1 = line([X0, X1], [Y0, Y1], [Z0, Z1], 'Color', "black", 'LineWidth',2);
        hold on;
        d2 = line([X1, X2], [Y1, Y2], [Z1, Z2], 'Color', "blue", 'LineWidth',2);
        hold on;
        d3 = line([X2, X3], [Y2, Y3], [Z2, Z3], 'Color', "yellow", 'LineWidth',2);
        hold on;
        d4 = line([X3, X4], [Y3, Y4], [Z3, Z4], 'Color', "red", 'LineWidth',2);
        hold on;
        axis([-600 400 -400 600 0 300]) 
        pause(0.05);
        delete(d1); delete(d2); delete(d3); delete(d4);
    end
    for t1=t1tr:-0.05:t1s
        t2 = t2s;
        l1 = l4+30;

        X1 = X0; Y1 = Y0; Z1 = Z0 + l1;
        X2 = X0 - sin(t1)*l2; Y2 = Y0 + cos(t1)*l2; Z2 = Z1;
        X3 = X2-sin(t1+t2)*l3; Y3 = Y2 + cos(t1+t2)*l3 ; Z3 = Z1;
        X4 = X3; Y4 = Y3; Z4 = Z3-l4;

        grid on
        d1 = line([X0, X1], [Y0, Y1], [Z0, Z1], 'Color', "black", 'LineWidth',2);
        hold on;
        d2 = line([X1, X2], [Y1, Y2], [Z1, Z2], 'Color', "blue", 'LineWidth',2);
        hold on;
        d3 = line([X2, X3], [Y2, Y3], [Z2, Z3], 'Color', "yellow", 'LineWidth',2);
        hold on;
        d4 = line([X3, X4], [Y3, Y4], [Z3, Z4], 'Color', "red", 'LineWidth',2);
        hold on;
        axis([-600 400 -400 600 0 300]) 
        pause(0.05);
        delete(d1); delete(d2); delete(d3); delete(d4);
    end

end
