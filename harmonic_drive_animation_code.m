clc
clear all  %clear all the variables


% Load the CSV file
profile = readmatrix("E:\temporary shift\SLP\tooth profile\tooth_profile.csv");   % change the file location accordingly

xcr = profile(1:30,1);  % X-coordinates of the right side of the tooth profile of circular spline
ycr = profile(1:30,2);  % Y-coordinates of the right side of the tooth profile of circular spline

xfr = profile(1:30,3);  % X-coordinates of the right side of the tooth profile of flexspline
yfr = profile(1:30,4);  % Y-coordinates of the right side of the tooth profile of flexspline

xcl = profile(31:60,1);  % X-coordinates of the left side of the tooth profile of circular spline
ycl = profile(31:60,2);  % Y-coordinates of the left side of the tooth profile of circular spline

xfl = profile(31:60,3);  % X-coordinates of the left side of the tooth profile of flexspline
yfl = profile(31:60,4);  % Y-coordinates of the left side of the tooth profile of flexspline


% parameters from the parameters file

%-----------------------------------------------------------------


% R_f0 = 77.151079136690640;
% w_0 = 1.653237410071943/2;
% r_0 = R_f0 - 2*w_0;       % prime circle radius
% %w_0 = (R_f0 - r_0)/2; % radial eccentricity 
% N_f = 280;           % number of teeth on flexspline
% N_c = N_f + 2;       % number of teeth on circular spline
% D_f0 = 153.2;        % diameter of undeformed flexspline

%--------------------------------------------------------------------

% parameters given below are from the 2011 paper, these gave close to accurate results
% NOTE: uncomment once the whole part between the lines to get the
% animation

%-------------------------------------------------------

r_0 = 60.31;       % prime circle radius
w_0 = 0.41713;       % accurate w_0 value till 5 decimal places
%w_0 = 0.42;         % value given in 2011 paper giving unaccurate results
%due to just 2 decimal places
N_f = 240;           % number of teeth on flexspline
N_c = N_f + 2;       % number of teeth on circular spline
D_f0 = 121.46;        % diameter of undeformed flexspline

rounds = 10;   % number of rounds of wave generator

% updating the profile for the above parameters
ycr = ycr - 75.4978 + 60.31;  
yfr = yfr - 75.4978 + 60.31;
ycl = ycl - 75.4978 + 60.31;  
yfl = yfl - 75.4978 + 60.31;

%----------------------------------------------------------------

% x_fM contains the x coordinates of a single tooth on the flexspline, both
% left and right included
% similarly y_fM for y coordinates on flexspline

x_fM(1:30) = flipud(xfr);   % reversing the order of coordinates for better organization
x_fM(31:60) = xfl;          % filling the xfl values in x_fM
y_fM(1:30) = flipud(yfr);   % similar for y coordinates
y_fM(31:60) = yfl;

R_M = [x_fM',y_fM'];        % combining the x and y coordinates of flexspline into a single 60x2 matrix

% x_CM contains the x coordinates of a single tooth on the flexspline, both
% left and right included
% similarly y_CM for y coordinates on flexspline

x_CM(1:30) = flipud(xcl);  % reversing the order of coordinates for better organization
x_CM(31:60) = xcr;         % filling the xfl values in x_fM
y_CM(1:30) = flipud(ycl);  % similar for y coordinates
y_CM(31:60) = ycr;

R_CM = [x_CM',y_CM'];      % combining the x and y coordinates of circular spline into a single 60x2 matrix
 
R = r_0 + w_0*(1+cos(2*0));  % eqaution for the cam in polar form solved at theta = 0
% R is the coordinate of the tooth base, point Of

R_of(1:60,1) = R*sin(0);    % x coordinate of the point Of (tooth base, as denoted in paper)
R_of(1:60,2) = R*cos(0);    % y coordinate of the point Of (tooth base, as denoted in paper)

r_M_of = R_M - R_of;        % r_M_of contains coordinates of point M (any point on the tooth profile) wrt to the point Of (tooth base)

theta_C = 0;                % angle of circular spline, kept zero beacuse 
% in the animation we fix circularspline and move wave genrator and observe
% the flexspline

% Goal: To find theta_f, given theta_W(input to wave generator)
% what happens in the lopp below? 
% 1)this is the loop interates over different values of wave generator, theta_W (input)
% 2)we use theta_W to find theta_CW, using theta_CW = theta_C - theta_W
% 3)We use this value of theta_CW to find theta (independent parameter) using equation 18, on page 3 of that paper
% 4)theta is then used to finally compute theta_f using relation: theta_f = theta_W + theta + phi;

theta_fss = zeros(180,1);


for j = 0:5:rounds*360    % range of angle of wave generator for which we are simulating
%for j = 119*360:5:120*360

    % STEP 1
    theta_W = deg2rad(j);      % from degree to radian
    
    % STEP 2
    theta_CW = theta_C - theta_W;    % relation between thses three variables
    
    % STEP 3
    syms theta_syms      % using symbolic variable for theta (independent variable)

    R_syms = r_0 + w_0*(1+cos(2*theta_syms));  % symbolic variable for R too
    R_dot_syms = 2*w_0*sin(2*theta_syms);      % derivative of R wrt theta

    f = sqrt(R_syms^2 + R_dot_syms^2); % defining some function f, which is the integrand in the relation of theta_CW and theta, written in just the next line
    %theta_CW = (2*N_f/(D_f0*N_c))*int(f,theta_syms,0,theta);
    % the relation of theta_CW and theta is written in the line above
    

    % Compute the integral symbolically, and denote it F_theta
    F_theta = matlabFunction(int(f, theta_syms, 0, theta_syms), 'Vars', theta_syms);

    % subtracting theta_CW from the integral and using fsolve to solve for
    % the theta where this "func_to_solve" becomes zero
    func_to_solve = @(theta) (2*N_f/(D_f0*N_c)) * F_theta(theta) - theta_CW;
    
    % Solve for theta using fzero
    theta_guess = 1;  % Initial guess
    theta = fzero(func_to_solve, theta_guess);

    R = r_0 + w_0*(1+cos(2*(theta)));   % equation of cam in polar form

    R_dot = -2*w_0*sin(2*(theta));      % derivative of R wrt to theta

    phi = acos(R/sqrt(R^2 + R_dot^2));  % phi as defined in paper, oscillating angle

    theta_f = theta_W + theta + phi;    % relation as in paper, page 4, eq 19

    R_of(1:60,1) = R*sin(theta_f - phi);    % updated x-coordinate of point Of (tooth base)
    R_of(1:60,2) = R*cos(theta_f - phi);    % updated y-coordinate of point Of (tooth base)
    % we subtract phi in the argument to align with the paper
    % refer to fig 7 on page 4 for better understanding

    B_fS = [cos(theta_f) sin(theta_f);            %rotation matrix from x,y to x',y'
           -sin(theta_f) cos(theta_f)];              

    % the above rotation matirc is to go from global coordinate system
    % S(O;X,Y) to S_f system
    
    % this for loop is for updating the tooth profile coordinates
    for i = 1:60
        R_M(i,:) = (R_of(i,:)' + B_fS*r_M_of(i,:)')';
        
    end

    % final animation plot for flexspline teeth
    plot(R_M(:,1),R_M(:,2),'r')
    axis equal
    hold on
    
    theta_fss(j/5+1) = theta;

    % the print statements below are just to check verify if the results
    % are correct, can be commented

    fprintf('theta_W = %.2f\n', rad2deg(theta_W));
    fprintf('theta = %.2f\n', rad2deg(theta));
    fprintf('theta_f = %.2f\n', rad2deg(theta_f));
    fprintf('phi = %.2f\n', rad2deg(phi));
    fprintf('error = %.2f\n', -rad2deg(theta_W)-rad2deg(theta)-360);
    disp("--------")

end



%% circular spline plot
for i = 0:2*rounds
    theta_C2 = deg2rad(360/N_c);       % angle at which second tooth of circular splien is present
    
     B = [cos(i*theta_C2) -sin(i*theta_C2);            %rotation matrix
          sin(i*theta_C2)  cos(i*theta_C2)];              
    
        % the above rotation matirx is to rotate the first tooth by  theta_C2
        % to get second tooth
    
     R_CM2 = zeros(60,2);     % storing the coordinates of the second tooth profile
    
    % this for loop is for updating the tooth profile coordinates
        for i = 1:60
            R_CM2(i,:) = B*R_CM(i,:)';
            
        end
    
    % plotting the circular spline tooth profiles
    %plot(R_CM(:,1), R_CM(:,2), 'b', 'LineWidth', 3)  
    hold on
    plot(R_CM2(:,1), R_CM2(:,2), 'b', 'LineWidth', 3)
    hold on
end

% CONCLUSION

% on using the parameter from the parameters metlab file, we get these
% values:

% theta_W = 180.00
% theta = -181.90
% theta_f = -1.82
% phi = 0.08

% the theta in this case was expected to be -181.2857, because in this
% reduction ratio was -1/140 (since N_f = 280), so in 360 degree rotation
% of wave generator, we should get theta(angle of point O_f) = 360/140 = 2.571428
% degrees more than theta_W in magnitude, so in our 180 case, 
% it should be -(180 + 2.571428/2) = -181.2857143 
% so we get an error of 0.6142857 degrees

% on using the parameter from the paper, we get these
% values:

% theta_W = 180.00
% theta = -181.48
% theta_f = -1.44
% phi = 0.04

% the theta in this case was expected to be -181.50, because in this
% reduction ratio was -1/120 (since N_f = 240), so in 360 degree rotation
% of wave generator, we should get theta(angle of point O_f) = 360/120 = 3
% degrees more than theta_W in magnitude, so in our 180 case, 
% it should be -(180 + 3/2) = -181.5 
% so we get an error of 0.02 degrees