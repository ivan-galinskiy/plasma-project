# A preliminary file to get started in writing the simulation

# The state of the system will be kept in two arrays. One will specify
# the charge/mass, positions and velocities of Np particles. 
# The positions are specified in 1D, while the velocities, for magnetic reasons,
# are specified in 3D.
# That is, it's an array with 6 columns and Np rows.

# Note that due to leapfrog integration, the particles' positions and 
# velocities ARE NOT SIMULTANEOUS. Specifically, the velocities "lag behind" by
# dt/2.

Np = 1;
Ng = 1000;
dx = 0.1;
L = 10;

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
    qom = particle_state[1];
    
    x = particle_state[2];
    
    vx = particle_state[3];
    vy = particle_state[4];
    vz = particle_state[5];
    
    a = field(world_grid, x)*qom;
    return [x vx+dt*a vy vz];
end

# A function to move, i.e. modify the position of a single particle
# (given by a row in "particles") according to its velocity.

# This function must be mapped with map! to the particles array

function move(particle_state, dt)
    qom = particle_state[1];
    
    x = particle_state[2];
    
    vx = particle_state[3];
    vy = particle_state[4];
    vz = particle_state[5];
    
    a = field(world_grid, x)*qom;
    return [x+dt*vx vx vy vz];
end

# A function to calculate the charge density on the grid. This is done by PIC.
# It uses the current particle positions from "particles"
