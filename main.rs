using PyPlot;
# A preliminary file to get started in writing the simulation

# The state of the system will be kept in two arrays. One will specify
# the charge/mass, positions and velocities of Np particles. 
# The positions are specified in 1D, while the velocities, for magnetic reasons,
# are specified in 3D.
# That is, it's an array with 6 columns and Np rows.

# Note that all the units are in CGS, such that the force between two charges
# (in 3D) is simply F = q1*q2/r^2. In 1D (charge sheets), this force is 
# F = q1*q2

# Note that due to leapfrog integration, the particles' positions and 
# velocities ARE NOT SIMULTANEOUS. Specifically, the velocities "lag behind" by
# dt/2.

Np = 1;
Ng = 10000;
dx = 0.1;
dt = 0.1;
L = 10.0;

particles = zeros(Np, 6);

# The other array world_grid will contain the densities and electric fields 
# specified at the grid points.
# The grid has Ng points

world_grid = zeros(Ng, 2);

# The initial conditions are specified as {x_i(0)}, {v_i(0)}. Note that this is
# inconsistent with the leapfrog scheme. Therefore, these will later be
# modified to fit in.

particles_initial = zeros(Np, 4);

# A function to accelerate, i.e. modify the velocity of a single particle
# (given by a row in "particles") according to the local field.

# This function must be mapped with map! to the particles array

function accelerate(particle_state, dt)
    q, m, x, vx, vy, vz = particle_state;
    
    a = field(world_grid, x)*qom;
    return [x vx+dt*a vy vz];
end

# A function to move, i.e. modify the position of a single particle
# (given by a row in "particles") according to its velocity.

# This function must be mapped with map! to the particles array

function move(particle_state, dt)
    q, m, x, vx, vy, vz = particle_state;
    
    xn = x+dt*vx;
    
    # Check if the new position if within bounds. If not, apply periodic
    # boundary conditions
    if xn > L
        xn = xn % L;
    elseif xn < 0
        xn = L - xn % L;
    end
    
    return [xn vx vy vz];
end

# A function to calculate the charge density on the grid. This is done by PIC.
# It uses the current particle positions from "particles".

# This function must be mapped with map to the particles array. It modifies
# the density values in "world_grid"

function rho_calculate(particle_state)
    q, m, x, vx, vy, vz = particle_state;
    
    # Point of the grid that is to the left of the particle
    left_point = Int(floor(x / dx)) + 1;
    
    # Consequently, the one to the right is
    right_point = Int(floor(x / dx)) + 2;
    
    # TODO: Add boundary conditions
    
    # The distance to the left point is
    dL = x % dx;
    
    # And to the right one
    dR = dx - dL;
    
    # Therefore, in the PIC scheme the density at the left point gains q*(dL/dx)
    world_grid[left_point, 1] += q*(dL/dx);
    
    # And analogously for the right point
    world_grid[right_point, 1] += q*(dR/dx);
    
    return;
end

function fields_calculate()
    # First we perform a Fast Fourier transform on the density
    rhok = fft(world_grid[:,1])
    
    # Then we obtain the potential for it (solving the Poisson equation)
    # TODO: fix the units and so on
    # TODO: implement the proper finite derivative
    phik = Array{Complex128}(length(rhok))
    phik[2:end] = [-rhok[k+1]/(k * (pi/L) * sinc(k/Ng))^2 for k in 1:length(rhok)-1]
    
    # And now we convert the potential back to x by IFFT
    world_grid[:,2] = real(ifft(phik))
    
    return;
end

# A temporary function to test the performance of the fields' calculation
function test_fields()
    # Let the density be a simple linear function
    for k in 1:length(world_grid[:,1])
        world_grid[k, 1] = 1;
    end
    
    # world_grid[5000, 1] = 1.0
    #world_grid[2000:2010, 1] = -1.0
    
    # Now we calculate the potential associated to it (which should be a 
    # cubic function)
    fields_calculate()
    
    # Let's see
    plot(world_grid[10:end-10,2])
end

test_fields()
