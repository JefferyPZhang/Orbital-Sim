M = 5.972e24; % mass of the Earth
earthRad = 6371; % Earth's radius (in km)
transparency = 0.2; % standard alpha value for surfaces
G = 6.67430e-11; % gravitational constant

% Display spherical coordinate system reference
sphRef = imshow("3D_Spherical.png");

% read file: picture of the earth (1024x512 dims)
original = imread('earth.png');
original = flipud(original);

% create sphere object
figure;
[x, y, z] = sphere(100);
surf(x * earthRad, y * earthRad, z * earthRad);

% mesh the picture onto the sphere
h = findobj('Type', 'surface');
set(h, 'CData', original, 'FaceColor', 'texture', 'edgecolor', 'none');

% visualization settings, labels, etc
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
set(gcf, 'Color', 'k');

% create dashed axis for the Earth
z_start = -12000;
z_end = 12000;
dash_length = 600;
num_dashes = (z_end - z_start) / dash_length;
z = linspace(z_start, z_end, num_dashes);
x = zeros(size(z));
y = zeros(size(z));
hold on;
plot3(x, y, z, 'b--', 'LineWidth', 2);
hold off;

% adjust view window
axis equal;
axis(earthRad * [-2 2 -2 2 -2 2]);

% save reference to the plot
myPlot = gcf;

valid = false
validMasses = [1, 2, 3, 6, 12];

% if not valid, ask again
while ~valid
    m = input("Enter a valid value for the mass (in kg, valid values include 1, 2, 3, 6, 12): ");
    in_vector = any(validMasses == m)
    if in_vector
        valid = true;
    end
end

valid = false

% if not valid, ask again
while ~valid
    r_0 = input("Enter a valid value for r (in km, between 150 & 2000 ABOVE the earth's surface): ");
    if isnumeric(r_0) && r_0 >= 150 && r_0 <= 2000
        valid = true;
    end
end
figure(myPlot);
hold on;

% add Earth's radius to get rho w.r.t. Earth's center
r_0 = r_0 + earthRad;

% create surface visualization
[x_sphere, y_sphere, z_sphere] = sphere(100);
sphereObject = surf(r_0 * x_sphere, r_0 * y_sphere, r_0 * z_sphere);
set(sphereObject, 'FaceAlpha', transparency, 'edgecolor', 'none');
hold off;

valid = false

% if not valid, ask again
while ~valid
    theta_0 = input("Enter a value for theta (in radians * pi / 6) - must be between 0 and 6: ");
    if isnumeric(theta_0) && theta_0 >= 0 && theta_0 <= 6
        valid = true;
    end
end
figure(myPlot);
hold on;

% poles are special cases, check if 
northPole = false;
southPole = false;
if (theta_0 == 0) 
    northPole = true;
    plot3(0, 0, r_0, 'ro', 'MarkerSize', 10);
    x_0 = 0;
    y_0 = 0;
    z_0 = r_0;
elseif (theta_0 == 6)
    theta_0 = 2 * pi;
    southPole = true;
    plot3(0, 0, -r_0, 'ro', 'MarkerSize', 10);
    x_0 = 0;
    y_0 = 0;
    z_0 = -r_0;
else
    % apply scaling value
    theta_0 = theta_0 * (pi / 3);

    % create surface visualization
    [r, phi] = meshgrid(linspace(0, r_0 * 2), linspace(0, 2 * pi));
    [X, Y, Z] = sph2cart(phi, theta_0, r);
    coneObject = surf(X, Y, Z, 'FaceAlpha', transparency, 'edgecolor', 'none');
end

valid = false
if (~northPole && ~southPole)
    % if not valid, ask again
    while ~valid
        phi_0 = input("Enter a value for phi (in radians * pi / 6) - must be between 0 and 12: ");
        if isnumeric(phi_0) && phi_0 >= 0 && phi_0 <= 12
            valid = true;
        end
    end
    figure(myPlot);
    hold on;

    % apply scaling value
    phi_0 = phi_0 * (pi / 3);

    % create surface visualization
    [X, Z] = meshgrid(linspace(0, 10000), linspace(-10000, 10000));
    Y = X * tan(phi_0);
    planeObject = surf(X, Y, Z, 'FaceAlpha', transparency, 'edgecolor', 'none');
    hold off;
end

if (~northPole && ~southPole)
    [x_0, y_0, z_0] = sph2cart(phi_0, theta_0, r_0);
    x_0 = -x_0;
    y_0 = -y_0;
    figure(myPlot);
    hold on;

    % plot the intersection point
    plot3(x_0, y_0, z_0, 'ro', 'MarkerSize', 10);
    hold off;
end

figure(myPlot);
hold on;
delete(coneObject);
delete(planeObject);

% create tangent plane to the sphere defined by rho at this point
plane_length = 10000;
plane_width = 10000;
x_min = x_0 - plane_length / 2;
x_max = x_0 + plane_length / 2;
y_min = y_0 - plane_width / 2;
y_max = y_0 + plane_width / 2;
[X, Y] = meshgrid(linspace(x_min, x_max), linspace(y_min, y_max));
Z = (z_0 ^ 2 - x_0 * (X - x_0) - y_0 * (Y - y_0)) / z_0;
secondPlaneObject = surf(X, Y, Z, 'FaceAlpha', transparency, 'EdgeColor', 'none');
set(sphereObject, 'FaceAlpha', transparency / 2);
N = [x_0, y_0, z_0]
normN = N / norm(N);

% create a directional vector that points toward the North pole while 
% still staying tangent to the sphere defined by rho
if (northPole || southPole)
    dir = [1, 0, 0]
else
    dir = [0, 0, earthRad] - N;
    dir_2 = (dot(dir, N) / (norm(N) ^ 2)) * N;
    dir = dir - dir_2;
    dir = dir / norm(dir);
end

% plot the directional vector
firstQuiver = quiver3(x_0, y_0, z_0, dir(1) * 5000, dir(2) * 5000, dir(3) * 5000, 'b', 'LineWidth', 2, ...
    'ShowArrowHead', 'on', 'MaxHeadSize', 20);
hold off;

figure(myPlot);
hold on;

% if not valid, ask again
valid = false;
while ~valid
    turnAngle = input("Enter a value for the turn angle (in radians * pi / 6) - must be between 0 and 12: ");
    if isnumeric(turnAngle) && turnAngle >= 0 && turnAngle <= 12
        valid = true;
    end
end
turnAngle = turnAngle * (pi / 6);
% find the rotated vector
if (northPole)
    dir = [cos(turnAngle), sin(turnAngle), 0];
elseif (southPole)
    dir = [cos(turnAngle), -sin(turnAngle), 0];
else
    w = cross(N, dir);
    dir = norm(dir) * ((cos(turnAngle) / norm(dir)) * dir + (sin(turnAngle) / norm(w)) * w);
    dir = dir / norm(dir);
end

% plot the rotated vector
delete(firstQuiver);
delete(secondPlaneObject);
rotatedQuiver = quiver3(x_0, y_0, z_0, dir(1) * 5000, dir(2) * 5000, dir(3) * 5000, 'b', 'LineWidth', 2, ...
    'ShowArrowHead', 'on', 'MaxHeadSize', 20);

% decrease the alpha of the sphere defined by rho, further
set(sphereObject, 'FaceAlpha', transparency / 5);
hold off;

% this is the velocity required to maintain a perfectly concentric orbit
% around the Earth
vels = sqrt((2 * G * M) / rho_0) * dir;

% fully defined state vector!
initialStateVector = [m; N'; vels';];

% ask user to start the simulation
valid = false;
while ~valid
    start = input("Enter 'start' to start the simulation: ", 's');
    if strcmpi(start, 'start')
        valid = true;
    end
end