%% Appendix 1
%% formation controller leader follower system
N = 7; %number of robots initialised 1 leader, 6 followers for each node of the hexagon
init_cond = generate_initial_conditions(N, 'Width', 2, 'Height', 2, 'Spacing', 0.5); %standard robotarium syntax for intialising 
%the robots in the simulator
RB = Robotarium('NumberOfRobots', N, 'ShowFigure', true, 'InitialConditions', init_cond);

iterations = 5000; %number of iterations that the simulation will run

fw =            [0	-1	 0	 2	 0	-1;...
                -1	 0	-1	 0	 2	 0;...
                 0	-1	 0	-1	 0 	 2;...
                 2 	 0	-1	 0	-1	 0;...  %follower Laplcian 6x6
                 0	 2	 0	-1 	 0	-1;...
                -1	 0	 2	 0	-1	 0];

L = zeros(N, N); %appending the leader to make a 7x7
L(2:N, 2:N) = fw;
L(2, 2) = L(2, 2) + 1;
L(2, 1) = -1;
disp(L); 

d = 0.4; %edge length for hexagon
ddiag = d*2; %diagonal length between agents
            
F_weights = [0      d	   0   ddiag     0       d;...
             d      0	   d     0     ddiag     0;...
             0      d	   0     d       0	  ddiag;... %hexaxgon distance weight definition wrt Laplacian
          ddiag 	0	   d     0       d	     0;...
             0	  ddiag	   0  	 d       0	     d;...
             d	    0	 ddiag	 0       d       0];
  
weights = zeros(N, N); 
weights(2:N, 2:N) = F_weights; % releasing leader from the hexagon formation 
disp(weights);

formation_control_gain = 10; %gain for the controller
dx = zeros(2, N); %iniialising follower velocities
state = 1; %initialising leader state to trael to waypoint 1

uni_barrier_cert = create_uni_barrier_certificate_with_boundary(); %calling a function for implementing the barrier certificates
si_to_uni_dyn = create_si_to_uni_dynamics('LinearVelocityGain', 0.8, 'AngularVelocityLimit', pi/2); %calling the unicylce dynamics converter
leader_controller = create_si_position_controller('XVelocityGain', 0.8, 'YVelocityGain', 0.8, 'VelocityMagnitudeLimit', 0.1); % initialising the robotarium single integrator position controller

way_pts = [-1 0.8; -1 -0.8; 1 -0.8; 1 0.8]'; %defining the waypoints for the leader to traverse through
close_enough = 0.05; % defining a waypoint tolerance


for t = 1:iterations
    
    x = RB.get_poses();% getting the initial positions of the agents from the simulator
    
    for i = 2:N %  initialising the control loop for the followers
        
        dx(:, i) = [0 ; 0];
        
        for j = topological_neighbors(L, i) %computing and updating the dx values wrt each neighbour
            dx(:, i) = dx(:, i) + formation_control_gain*(norm(x(1:2, i) - x(1:2, j))^2 - weights(i, j)^2)*(x(1:2, j) - x(1:2, i));
            %the controller
        end 
    end
    
    way_pt = way_pts(:, state);
    %running the leader through the waypoints
    switch state        
        case 1
            dx(:, 1) = leader_controller(x(1:2, 1), way_pt);
            if(norm(x(1:2, 1) - way_pt) < close_enough)
                state = 2;
            end
        case 2
            dx(:, 1) = leader_controller(x(1:2, 1), way_pt);
            if(norm(x(1:2, 1) - way_pt) < close_enough)
                state = 3;
            end
        case 3
            dx(:, 1) = leader_controller(x(1:2, 1), way_pt);
            if(norm(x(1:2, 1) - way_pt) < close_enough)
                state = 4;
            end
        case 4
            dx(:, 1) = leader_controller(x(1:2, 1), way_pt);
            if(norm(x(1:2, 1) - way_pt) < close_enough)
                state = 1;
            end
    end
    
    %standard robotarium syntax to limit actuator speed for protection
    norms = arrayfun(@(x) norm(dx(:, x)), 1:N);
    threshold = 3/4*RB.max_linear_velocity;
    to_thresh = norms > threshold;
    dx(:, to_thresh) = threshold*dx(:, to_thresh)./norms(to_thresh);
    
    %converting the single integrator commands into unicycle xyTheta
    %commands and running them through the barrier certificates for updated
    %collision free trajectory
    dxu = si_to_uni_dyn(dx, x);  
    dxu = uni_barrier_cert(dxu, x);
    
    %pushing velocities to the agents in the simulator
    RB.set_velocities(1:N, dxu);
    
    RB.step(); 
end

%standard simulator code for error checks
RB.debug();
