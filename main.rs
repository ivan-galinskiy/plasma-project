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

Np = Int(floor((2^6-1)*70)); # Number of particles
Ng = Int(2^6); # Number of grid points (power of 2 for efficiency)
Nt = Int(500); # Number of time points

L = 100.0; # Length of the simulation cell (in cm)
#L = 0.07*4*pi*Np; # Length of the simulation cell (in cm)
T = 100.0; # Time length of the simulation

dx = float(L/(Ng-1)); # Size of each grid interval
dt0 = float(T/Nt);

particles = zeros(Np, 6);

# The other array world_grid will contain the densities and electric fields 
# specified at the grid points.
# The grid has Ng points (for Ng-1 intervals), plus an additional point to
# make periodic boundaries work

world_grid = zeros(Ng+1, 3);

# The initial conditions are specified as {x_i(0)}, {v_i(0)}. Note that this is
# inconsistent with the leapfrog scheme. Therefore, these will later be
# modified to fit in.

particles_initial = zeros(Np, 6);

# A function to accelerate, i.e. modify the velocity of a single particle
# (given by a row in "particles") according to the local field.

# This function must be mapped with map! to the particles array

function null_world()
    particles = zeros(Np, 6);
    world_grid = zeros(Ng+1, 3);
    particles_initial = zeros(Np, 6);
end

function accelerate(dt=1e-3)
    for n in 1:Np
        q, m, x, vx, vy, vz = particles[n,:];
    
        a = field_interpolate(x)*q/m ;
        particles[n,:] = [q m x vx+dt*a vy vz];
    end
end

function field_calculate()
    for n in 2:Ng
        world_grid[n, 3] = (world_grid[n-1, 2] - world_grid[n+1, 2])/(2*dx)
    end
    world_grid[1, 3] = (world_grid[Ng, 2] - world_grid[2, 2])/(2*dx)
    world_grid[Ng+1, 3] = world_grid[1, 3]
end

function field_interpolate(x)
    N = x/dx + 1
    left_point = Int(floor(N))
    right_point = left_point + 1
    
    dL = x - dx*(left_point - 1)
    
    return world_grid[left_point, 3] + (world_grid[right_point, 3] - world_grid[left_point, 3])*dL/dx
end

# A function to move, i.e. modify the position of a single particle
# (given by a row in "particles") according to its velocity.

# This function must be mapped with map! to the particles array


function move(dt=1e-3)
    for n in 1:Np
        q, m, x, vx, vy, vz = particles[n,:];
        
        x += dt*vx;
        
        # Check if the new position is within bounds. If not, apply periodic
        # boundary conditions
        if x >= L
            x -= L;
        elseif x < 0
            x += L;
        end
        
        particles[n,:] = [q m x vx vy vz];
    end
end

# A function to calculate the charge density on the grid. This is done by PIC.
# It uses the current particle positions from "particles".

# It modifies the density values in "world_grid"

function rho_calculate()
    world_grid[:,1] = 0
    for n in 1:Np
        q, m, x, vx, vy, vz = particles[n,:];
        
        if x >= L
            x -= L;
        end
        
        # Point of the grid that is to the left of the particle
        left_point = Int(floor(x / dx))+1;
        
        # Consequently, the one to the right is
        right_point = left_point + 1;
        
        # The distance to the left point is
        dL = x - (left_point-1) * dx;
        
        # And to the right one
        dR = dx - dL;
        
        # Therefore, in the PIC scheme the density at the left point gains q*(dL/dx)
        world_grid[left_point, 1] += q*(1 - dL/dx)/dx;
        
        # And analogously for the right point
        world_grid[right_point, 1] += q*(1 - dR/dx)/dx;
    end
    
    # And we apply the periodic boundary conditions
    ld = world_grid[1, 1]
    rd = world_grid[Ng, 1]
    world_grid[1, 1] = ld + rd
    world_grid[Ng, 1] = world_grid[1, 1]
    
    return;
end

# For periodic boundaries, the total charge MUST be 0 in order for Maxwell
# equations to hold.
function potential_calculate()
    # First we perform a Fast Fourier transform on the density
    # Note that we omit the last point because it's periodic
    
    world_grid[1, 1] += world_grid[Ng+1, 1]
    world_grid[Ng+1, 1] = world_grid[1, 1]
    
    rhok = fft(world_grid[1:Ng,1])
    
    # Then we obtain the potential for it (solving the Poisson equation).
    # The 4pi factor is due to our use of the CGS system
    phik = -4*pi*rhok
    
    # This should be zero, so it was probably wise to remove this.
    #a0 = phik[1]
    phik[1] = 0
    
    function K(k)
        return 1/(k * (2*pi/L) * sinc(k/Ng))^2
    end
    
    for k in 2:Int(Ng/2)
        phik[k] = -phik[k] * K(k-1)
        phik[Ng + 2 - k] = -phik[Ng + 2 - k] * K(k-1)
    end
    phik[Ng/2+1] = -phik[Ng/2+1] * K(Ng/2)
    
    # And now we convert the potential back to x by IFFT
    world_grid[1:Ng,2] = real(ifft(phik))
    
    # And apply boundary conditions
    world_grid[Ng+1,2] = world_grid[1,2]
    
    return;
end

# Populate the simulation space uniformly with two streams of same species with
# opposing velocities (+- v). Np assumed to be odd.
function two_streams_populate(q, m, v)
    for n in 1:Np
        particles[n,:] = [q m (L * (n-1)/Np ) v*(-1)^n 0 0];
    end
end

# A temporary function to test the performance of the fields' calculation
function test_fields()
    world_grid[Int(floor(Ng/5)):Int(floor(Ng/5))+1, 1] =  1.0
    world_grid[Int(floor(Ng/4)):Int(floor(Ng/4))+1, 1] = -1.0
    
    # Now we calculate the potential associated to it.
    potential_calculate()
    
    # Let's see
    figure()
    plot(linspace(0, L, Ng), world_grid[1:Ng,2])
    title("Potential of two opposite delta charges")
    
    figure()
    field_calculate()
    plot(linspace(0, L, Ng), world_grid[1:Ng,3])
    title("Field of two opposite delta charges")
    #plot(xr, [field_interpolate(x) for x in xr])
end

function test_rho()
    particles[1,:] = [1 1 L/5 0 0 0]
    particles[2,:] = [-1 1 L/4 0 0 0]
    rho_calculate()
    
    figure()
    plot(linspace(0, L, Ng), world_grid[1:Ng,1])
    title("Density of two opposite point charges")
    
    null_world()
    two_streams_populate(1, 1, 0)
    
    rho_calculate()
    potential_calculate()
    
    figure()
    plot(linspace(0, L, Ng), world_grid[1:Ng,1], "o")
    title("Density of a uniform charge distribution")
end

function test_loop()
    particles[1,:] = [1 0.1 L/2-L/10 0 0 0];
    particles[2,:] = [-1 0.1 L/4+L/10 0 0 0];
    
    # All the other particles are not important
    for k in 3:Np
        particles[k,:] = [0 1 0 0 0 0]
    end
    
    rho_calculate()
    potential_calculate()
    field_calculate()
    
    #plot(world_grid[1:Ng,3])
    #println(err[end])
    #yscale("log")
    xr = linspace(0, L, 10000)
    #plot(xr, [field_interpolate(x) for x in xr])
    #plot(world_grid[:,1])
    #plot(world_grid[:,2])
    
    # Advance the velocities half a step into the past for leapfrog
    accelerate(-dt0/2)
    
    pos1 = zeros(Nt)
    pos2 = zeros(Nt)
    
    for t in 1:Nt
        rho_calculate()
        potential_calculate()
        field_calculate()
        
        accelerate(dt0)
        move(dt0) 
        pos1[t] = particles[1,3]
        pos2[t] = particles[2,3]
    end
    
    pos_dif = pos1 - pos2
    plot(pos_dif)
    
    # Calculate the second derivative in time to obtain the force
    pos_der = zeros(length(pos_dif))
    for k in 2:length(pos_dif)-1
        pos_der[k] = (pos_dif[k-1] + pos_dif[k+1] - 2*pos_dif[k])/dt0^2
    end
    pos_der[1] = pos_der[2]
    pos_der[end] = pos_der[end-1]
    
    figure()
    plot(pos_der)
end

# Density history
dens_his = zeros(Nt, Ng)

# Phase space history
phase_his = zeros(Nt, 2*Np)

function two_streams_test()
    null_world()
    two_streams_populate(1, 0.1, 1)
    
    # Now we induce a perturbation in the particles' velocities
    for n in 1:Np
        q, m, x, vx, vy, vz = particles[n,:];
        vxn = 0;
        if vx > 0
            vxn = vx + 0.2*sin(2*pi*x/L)
        else
            vxn = vx - 0.2*sin(2*pi*x/L)
        end
        
        particles[n,:] = [q, m, x, vxn, vy, vz]
    end
    
    # Advance the velocities half a step into the past for leapfrog
    rho_calculate()
    potential_calculate()
    field_calculate()
    accelerate(-dt0/2)
    
    for t in 1:Nt
        rho_calculate()
        potential_calculate()
        field_calculate()
        
        dens_his[t,:] = world_grid[1:Ng, 1]
        phase_his[t, 1:Np] = particles[:,3]
        phase_his[t, Np+1:2*Np] = particles[:,4]
        
        accelerate(dt0)
        move(dt0)
    end
end

# Force compilation
two_streams_test()
