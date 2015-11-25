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
Ng = 128; # Number of grid points

L = 1.0;

dx = L/(Ng-1);
dt0 = 1e-4;

particles = zeros(Np, 6);

# The other array world_grid will contain the densities and electric fields 
# specified at the grid points.
# The grid has Ng points (for Ng-1 intervals)

world_grid = zeros(Ng, 3);

# The initial conditions are specified as {x_i(0)}, {v_i(0)}. Note that this is
# inconsistent with the leapfrog scheme. Therefore, these will later be
# modified to fit in.

particles_initial = zeros(Np, 6);

# A function to accelerate, i.e. modify the velocity of a single particle
# (given by a row in "particles") according to the local field.

# This function must be mapped with map! to the particles array

function accelerate(dt=1e-3)
    for n in 1:Np
        q, m, x, vx, vy, vz = particles[n,:];
    
        a = field_interpolate(x)*q/m ;
        #println(a)
        particles[n,:] = [q m x vx+dt*a vy vz];
    end
end

function field_calculate()
    world_grid[1, 3] = (world_grid[end-1, 2] - world_grid[2, 2])/(2*dx)
    world_grid[end, 3] = world_grid[1, 3]
    
    for n in 2:Ng-1
        world_grid[n, 3] = (world_grid[n-1, 2] - world_grid[n+1, 2])/(2*dx)
    end
end

function field_interpolate(x)
    N = x/dx + 1
    left_point = Int(floor(N))
    right_point = left_point + 1
    if left_point == Ng
        right_point = 2
    end
    
    dL = x - dx*(left_point - 1)
    
    return world_grid[left_point, 3] + (world_grid[right_point, 3] - world_grid[left_point, 3])*dL/dx
end

# A function to move, i.e. modify the position of a single particle
# (given by a row in "particles") according to its velocity.

# This function must be mapped with map! to the particles array


function move(dt=1e-3)
    for n in 1:Np
        q, m, x, vx, vy, vz = particles[n,:];
        
        xn = x+dt*vx;
        
        # Check if the new position if within bounds. If not, apply periodic
        # boundary conditions
        if xn >= L
            xn = xn - L;
        elseif xn < 0
            xn = L + xn;
        end
        
        particles[n,:] = [q m xn vx vy vz];
    end
end

# A function to calculate the charge density on the grid. This is done by PIC.
# It uses the current particle positions from "particles".

# It modifies the density values in "world_grid"

function rho_calculate()
    world_grid[:,1] = 0
    for n in 1:Np
        q, m, x, vx, vy, vz = particles[n,:];
        
        # Point of the grid that is to the left of the particle
        left_point = Int(floor(x / dx))+1;

        if left_point == Ng
            left_point = Ng-1
        end
        
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
    world_grid[1, 1] += world_grid[end, 1]
    world_grid[end, 1] = world_grid[1, 1]
    
    return;
end

# For periodic boundaries, the total charge MUST be 0 in order for Maxwell
# equations to hold.
function potential_calculate()
    # First we perform a Fast Fourier transform on the density
    rhok = fft(world_grid[:,1])
    
    # Then we obtain the potential for it (solving the Poisson equation).
    # The 4pi factor is due to our use of the CGS system
    phik = -4*pi*rhok
    
    # This should be zero, so it was probably wise to remove this.
    #a0 = phik[1]
    phik[1] = 0
    
    for k in 1:(Ng-1)
        #phik[k+1] = -phik[k+1] / (k * (2*pi/L) * sinc(k/Ng))^2 + 0.5*a0*(im*L/(2*pi*k))^2
        phik[k+1] = -phik[k+1] / (k * (2*pi/L) * sinc(k/Ng))^2
    end
    
    # And now we convert the potential back to x by IFFT
    world_grid[:,2] = real(ifft(phik))
    
    return;
end

# A temporary function to test the performance of the fields' calculation
function test_fields()
    world_grid[Int(floor(Ng/3)):Int(floor(Ng/3))+1, 1] =  1.0
    world_grid[Int(floor(2*Ng/3)):Int(floor(2*Ng/3))+1, 1] = -1.0
    
    # Now we calculate the potential associated to it.
    potential_calculate()
    
    # Let's see
    #plot(world_grid[1:end,2])
    field_calculate()
    #plot(world_grid[10:end-10,3])
    xr = linspace(0, L, 10000)
    plot(xr, [field_interpolate(x) for x in xr])
end

function test_rho()
    #for n in 1:Np
    #    particles[n,:] = [1 1 0.9999*n*L/Np 0 0 0]
    #end
    particles[1,:] = [1 1 L/2 0 0 0]
    particles[2,:] = [-1 1 L/2+1e-6 0 0 0]
    rho_calculate()
    plot(world_grid[:,1])
end

function test_loop()
    
    err = zeros(1000)
    for k in 1:1000
        px = k*L/1000
        particles[1,:] = [1 1 px 0.1 0 0];
        #particles[2,:] = [1 1 L/2+1 0 0 0];
        
        rho_calculate()
        # Add a neutralizing background
        #temp = fft(world_grid[:,1])
        #temp[1,1] = 0
        #world_grid[:,1] = real(ifft(temp))
        
        #world_grid[:, 1] -= 2/L
        potential_calculate()
        field_calculate()
        err[k] = field_interpolate(px)
    end
    plot(err[:])
    yscale("log")
    xr = linspace(0, L, 10000)
    #plot(xr, [field_interpolate(x) for x in xr])
    #plot(world_grid[:,1])
    #plot(world_grid[20:40,3], "o")
    
%;#    accelerate(-dt0/2)
%;#    
%;#    Nt = 200
%;#    pos1 = zeros(Nt)
%;#    pos2 = zeros(Nt)
%;#    
%;#    for t in 1:Nt
%;#        rho_calculate()
%;#        potential_calculate()
%;#        field_calculate()
%;#        
%;#        
%;#        pos1[t] = field_interpolate(particles[1,3])
%;#        
%;#        accelerate(dt0)
%;#        move(dt0) #particles[1,3]
%;#        #pos2[t] = particles[2,3]
%;#    end
%;#    
%;#    plot(pos1[1:end])
%;#    #plot(pos2[1:end])
    
end

#test_fields()
test_loop()
#test_rho()
