%% Appendix 2
%% Hexagon formation control

N = 6;
init_cond = generate_initial_conditions(N, 'Width', 2, 'Height', 2, 'Spacing', 0.5);
RB = Robotarium('NumberOfRobots', N, 'ShowFigure', true, 'InitialConditions', init_cond);

formation_control_gain = 10;

iterations = 5000;

L =     [0	-1	 0	 2	 0	-1;...
        -1	 0	-1	 0	 2	 0;...
         0	-1	 0	-1	 0 	 2;...
         2 	 0	-1	 0	-1	 0;...
         0	 2	 0	-1 	 0	-1;...
        -1	 0	 2	 0	-1	 0];
 
d = 0.4; 

ddiag = d*2; 
            
weights = [0      d	   0   ddiag     0       d;...
             d      0	   d     0     ddiag     0;...
             0      d	   0     d       0	  ddiag;...
          ddiag 	0	   d     0       d	     0;...
             0	  ddiag	   0  	 d       0	     d;...
             d	    0	 ddiag	 0       d       0];
  

     
dx = zeros(2, N);
state = 1;
 

uni_barrier_cert = create_uni_barrier_certificate_with_boundary();
si_to_uni_dyn = create_si_to_uni_dynamics('LinearVelocityGain', 0.8, 'AngularVelocityLimit', pi/2);
leader_controller = create_si_position_controller('XVelocityGain', 0.8, 'YVelocityGain', 0.8, 'VelocityMagnitudeLimit', 0.1);

waypoints = [-0.75 0.75; -0.2 0.7; 0.8 0; 0.2 -0.4]';
close_enough = 0.05;

oone = 1;
for t = 1:iterations
    
    x = RB.get_poses();
    
    for i = 1:N
        
        dx(:, i) = [0 ; 0];
        
        for j = topological_neighbors(L, i)
            

            dx(:, i) = dx(:, i) + formation_control_gain*(norm(x(1:2, i) - x(1:2, j))^2 - weights(i, j)^2)*(x(1:2, j) - x(1:2, i));
            
        
            
        end 
    end
    
    norms = arrayfun(@(x) norm(dx(:, x)), 1:N);
    threshold = 3/4*RB.max_linear_velocity;
    to_thresh = norms > threshold;
    dx(:, to_thresh) = threshold*dx(:, to_thresh)./norms(to_thresh);
    
    dxu = si_to_uni_dyn(dx, x);  
    dxu = uni_barrier_cert(dxu, x);
    
    RB.set_velocities(1:N, dxu);
    
    RB.step(); 
    
end

r.debug();
