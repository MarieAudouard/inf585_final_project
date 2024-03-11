#include "simulation.hpp"

using namespace cgp;

// Convert a density value to a pressure
float density_to_pressure(float rho, float rho0, float stiffness)
{
	return stiffness*(rho-rho0);
}

float W_laplacian_viscosity(vec3 const& p_i, vec3 const& p_j, float h)
{
    // To do ...
    //  Fill it with laplacian of W_viscosity
    float const r = norm(p_i - p_j);
    assert_cgp_no_msg(r <= h);
    return 45.0 / (3.14159f * std::pow(h, 6)) * (h - r);
    //return 0.0f;
}

vec3 W_gradient_pressure(vec3 const& p_i, vec3 const& p_j, float h)
{
    // To do ...
    //  Fill it with gradient of W_spiky
    float const r = norm(p_i - p_j);
    assert_cgp_no_msg(r <= h);
    return - 45.0f / (3.14159f * std::pow(h, 6)) * std::pow(h - r, 2.0f) / r * (p_i-p_j);
    //return (p_i-p_j)/norm(p_i-p_j);
}

float W_density(vec3 const& p_i, const vec3& p_j, float h)
{
	float const r = norm(p_i-p_j);
    assert_cgp_no_msg(r<=h);
	return 315.0/(64.0*3.14159f*std::pow(h,9)) * std::pow(h*h-r*r, 3.0f);
}


void update_density(numarray<particle_element>& particles, float h, float m)
{
    // To do: Compute the density value (particles[i].rho) at each particle position
    //  rho_i = \sum_j m W_density(pi,pj)
    int const N = particles.size();
    for (int i = 0; i < N; ++i) {
        //particles[i].rho = 1.0f; // to be modified
        particles[i].rho = 0.0f;
        for (int j = 0; j < N; ++j)
        {
            if (norm(particles[i].p - particles[j].p) < h) {
                particles[i].rho += m * W_density(particles[i].p, particles[j].p, h);
            }
        }
    }
}

// Convert the particle density to pressure
void update_pressure(numarray<particle_element>& particles, float rho0, float stiffness)
{
	const int N = particles.size();
    for(int i=0; i<N; ++i)
        particles[i].pressure = density_to_pressure(particles[i].rho, rho0, stiffness);
}

// Compute the forces and update the acceleration of the particles
void update_force(numarray<particle_element>& particles, float h, float m, float nu)
{
	// gravity
    const int N = particles.size();
    for (int i = 0; i < N; ++i) {
        particles[i].f = m * vec3{ 0,-9.81f,0 };

        //TO Do
        // For all particles i
        //   Compute F_pressure
        //   Compute F_viscosity
        //   particles[i].f += (F_pressure + F_viscosity)
        // ...

        vec3 F_pressure = { 0,0,0 };
        for (int j = 0; j < N; ++j) {
			if (i != j && norm(particles[i].p - particles[j].p) < h) {
				F_pressure += -m * m * (particles[i].pressure + particles[j].pressure) / (2 * particles[j].rho * particles[i].rho) * W_gradient_pressure(particles[i].p, particles[j].p, h);
			}
		}
		particles[i].f += F_pressure;

		vec3 F_viscosity = { 0,0,0 };
        for (int j = 0; j < N; ++j) {
            if (i != j && norm(particles[i].p - particles[j].p) < h) {
                F_viscosity += nu * m * m * (particles[j].v - particles[i].v) / particles[j].rho * W_laplacian_viscosity(particles[i].p, particles[j].p, h);
            }
        }
        particles[i].f += F_viscosity;

    }
}

void simulate(float dt, numarray<particle_element>& particles, sph_parameters_structure const& sph_parameters)
{

	// Update values
    update_density(particles, sph_parameters.h, sph_parameters.m);                   // First compute updated density
    update_pressure(particles, sph_parameters.rho0, sph_parameters.stiffness);       // Compute associated pressure
    update_force(particles, sph_parameters.h, sph_parameters.m, sph_parameters.nu);  // Update forces

	// Numerical integration
	float const damping = 0.005f;
    int const N = particles.size();
	float const m = sph_parameters.m;
	for(int k=0; k<N; ++k)
	{
		vec3& p = particles[k].p;
		vec3& v = particles[k].v;
		vec3& f = particles[k].f;

		v = (1-damping)*v + dt*f/m;
		p = p + dt*v;
	}


	// Collision
    float const epsilon = 1e-3f;
    for(int k=0; k<N; ++k)
    {
        vec3& p = particles[k].p;
        vec3& v = particles[k].v;

        // small perturbation to avoid alignment
        if( p.y<-1 ) {p.y = -1+epsilon*rand_uniform();  v.y *= -0.5f;}
        if( p.x<-1 ) {p.x = -1+epsilon*rand_uniform();  v.x *= -0.5f;}
        if( p.x>1 )  {p.x =  1-epsilon*rand_uniform();  v.x *= -0.5f;}
    }

}

void solve_incompressibility_cell(grid_structure& grid, int i, int j, int k) {
    int N = grid.cellnumber;
    int nb_neighbours = 0;
    float div = 0;
    
    if (i > 0) {
        // to check that it is not in the ground
        div -= grid.velocity_x(i, j, k);
        nb_neighbours++;
    }
    if (i < N - 1) {
        div += grid.velocity_x(i + 1, j, k);
        nb_neighbours++;
	}

    if (j > 0) {
        // to check that it is not in the ground
        div -= grid.velocity_y(i, j, k);
        nb_neighbours++;
    }
    if (j < N - 1) {
		div += grid.velocity_y(i, j + 1, k);
		nb_neighbours++;
	}

    if (k > 0) {
        // to check that it is not in the ground
        div -= grid.velocity_z(i, j, k);
        nb_neighbours++;
    }
    if (k < N - 1) {
        div += grid.velocity_z(i, j, k + 1);
        nb_neighbours++;
    }

    if (nb_neighbours > 0 && div != 0) {
		div /= nb_neighbours;

        // we keep velocity = 0 for the ground/walls
        if (i > 0) {
            grid.velocity_x(i, j, k) += div;
        }
        if (j > 0) {
            grid.velocity_y(i, j, k) += div;
        }
        if (k > 0) {
            grid.velocity_z(i, j, k) += div;
        }
        if (i < N - 1) {
            grid.velocity_x(i + 1, j, k) -= div;
        }
        if (j < N - 1) {
			grid.velocity_y(i, j + 1, k) -= div;
		}
        if (k < N - 1) {
			grid.velocity_z(i, j, k + 1) -= div;
		}
	}
}

void solve_incompressibility(grid_structure& grid) {
    int nb_iter = 10; // gauss-seidel iterations
    for (int iter = 0; iter < nb_iter; ++iter) {
        for (int i = 0; i < grid.cellnumber; ++i) {
            for (int j = 0; j < grid.cellnumber; ++j) {
                for (int k = 0; k < grid.cellnumber; ++k) {
                    solve_incompressibility_cell(grid, i, j, k);
				}
			}
        }
    }
}

void update_grid(cgp::numarray<particle_element> const& particles, grid_structure& grid, pic_parameters_structure const& sph_parameters, float dt) {
    // set zero everywhere
    grid.velocity_x.fill(0);
    grid.velocity_y.fill(0);
    grid.velocity_z.fill(0);

    grid.empty.fill(emptycell);

    cgp::grid_3D<float> weights_x(grid.cellnumber);
    cgp::grid_3D<float> weights_y(grid.cellnumber);
    cgp::grid_3D<float> weights_z(grid.cellnumber);
    weights_x.fill(0);
    weights_y.fill(0);
    weights_z.fill(0);

    // get grid speed from particles
    for (int i = 0; i < particles.size(); i++) {
        particle_element const& part = particles[i];
        int coords_bottom_left[3] = { floor(part.p[0] / grid.cellsize), floor(part.p[1] / grid.cellsize), floor(part.p[2] / grid.cellsize) };
        float const dx = (part.p[0] / grid.cellsize) - coords_bottom_left[0];
        float const dy = (part.p[1] / grid.cellsize) - coords_bottom_left[1];
        float const dz = (part.p[2] / grid.cellsize) - coords_bottom_left[2];

        // update velocities from particles
        // we don't update when the grid is in the ground
        if (coords_bottom_left[0] > 0) {
            grid.velocity_x(coords_bottom_left[0], coords_bottom_left[1], coords_bottom_left[2]) += dx * part.v[0];
            weights_x(coords_bottom_left[0], coords_bottom_left[1], coords_bottom_left[2]) += dx;
            grid.empty(coords_bottom_left[0], coords_bottom_left[1], coords_bottom_left[2]) = fullcell;
        }
        if (coords_bottom_left[1] > 0) {
			grid.velocity_y(coords_bottom_left[0], coords_bottom_left[1], coords_bottom_left[2]) += dy * part.v[1];
			weights_y(coords_bottom_left[0], coords_bottom_left[1], coords_bottom_left[2]) += dy;
		}
        if (coords_bottom_left[2] > 0) {
			grid.velocity_z(coords_bottom_left[0], coords_bottom_left[1], coords_bottom_left[2]) += dz * part.v[2];
			weights_z(coords_bottom_left[0], coords_bottom_left[1], coords_bottom_left[2]) += dz;
		}

        if (coords_bottom_left[0] < grid.cellnumber - 1) {
            grid.velocity_x(coords_bottom_left[0] + 1, coords_bottom_left[1], coords_bottom_left[2]) += (1-dx) * part.v[0];
            weights_x(coords_bottom_left[0] + 1, coords_bottom_left[1], coords_bottom_left[2]) += (1-dx);
        }
        if (coords_bottom_left[1] < grid.cellnumber - 1) {
			grid.velocity_y(coords_bottom_left[0], coords_bottom_left[1] + 1, coords_bottom_left[2]) += (1-dy) * part.v[1];
            weights_y(coords_bottom_left[0], coords_bottom_left[1] + 1, coords_bottom_left[2]) += (1-dy);
		}
        if (coords_bottom_left[2] < grid.cellnumber - 1) {
            grid.velocity_z(coords_bottom_left[0], coords_bottom_left[1], coords_bottom_left[2] + 1) += (1-dz) * part.v[2];
            weights_z(coords_bottom_left[0], coords_bottom_left[1], coords_bottom_left[2] + 1) += (1-dz);
        }
    }


    // update velocities with gravity and use weights to normalize
    for (int i = 0; i < grid.cellnumber; ++i) {
		for (int j = 0; j < grid.cellnumber; ++j) {
			for (int k = 0; k < grid.cellnumber; ++k) {
                if (weights_x(i, j, k) > 0) {
                    grid.velocity_x(i, j, k) /= weights_x(i, j, k);
                }
                if (weights_y(i, j, k) > 0) {
					grid.velocity_y(i, j, k) /= weights_y(i, j, k);
				}
                if (weights_z(i, j, k) > 0) {
					grid.velocity_z(i, j, k) /= weights_z(i, j, k);
				}
                if (j > 0) {
					grid.velocity_y(i, j, k) += -9.81f * dt; // gravity, we keep velocity = 0 for the ground
                }
			}
		}
	}
    solve_incompressibility(grid);
}

void update_particle_velocities(cgp::numarray<particle_element>& particles, grid_structure const& grid, pic_parameters_structure const& sph_parameters, float dt) {
	// update velocities
	for (int i = 0; i < particles.size(); ++i) {
		// update velocity
        vec3 p = particles[i].p;
        vec3 v = particles[i].v;
		v = vec3(0,0,0);

        int coords_bottom_left[3] = { floor(p[0] / grid.cellsize), floor(p[1] / grid.cellsize), floor(p[2] / grid.cellsize) };
        float bar_coords[3] = { (p[0] / grid.cellsize) - coords_bottom_left[0], (p[1] / grid.cellsize) - coords_bottom_left[1], (p[2] / grid.cellsize) - coords_bottom_left[2]};

        v[0] = bar_coords[0] * grid.velocity_x(coords_bottom_left[0], coords_bottom_left[1], coords_bottom_left[2]);
        v[1] = bar_coords[1] * grid.velocity_y(coords_bottom_left[0], coords_bottom_left[1], coords_bottom_left[2]);
        v[2] = bar_coords[2] * grid.velocity_z(coords_bottom_left[0], coords_bottom_left[1], coords_bottom_left[2]);
        if (coords_bottom_left[0] < grid.cellnumber - 1) {
			v[0] += (1 - bar_coords[0]) * grid.velocity_x(coords_bottom_left[0] + 1, coords_bottom_left[1], coords_bottom_left[2]);
		}
        if (coords_bottom_left[1] < grid.cellnumber - 1) {
            v[1] += (1 - bar_coords[1]) * grid.velocity_y(coords_bottom_left[0], coords_bottom_left[1] + 1, coords_bottom_left[2]);
        }
        if (coords_bottom_left[2] < grid.cellnumber - 1) {
			v[2] += (1 - bar_coords[2]) * grid.velocity_z(coords_bottom_left[0], coords_bottom_left[1], coords_bottom_left[2] + 1);
		}

		particles[i].v = v;
	}
}

void update_particle_positions(cgp::numarray<particle_element>& particles, grid_structure& grid, pic_parameters_structure const& pic_parameters, float dt) {
	// set empty everywhere
    grid.empty.fill(emptycell);
    float const epsilon = 1e-3f;

    // update positions
	for (int i = 0; i < particles.size(); ++i) {
		vec3 p = particles[i].p;
		vec3 v = particles[i].v;
		vec3 f = particles[i].f;

		// update position
		p += dt * v;

        // Collision
        // small perturbation to avoid alignment
        if (p.x < 0) { p.x = epsilon * rand_uniform();  v.x *= -0.5f; }
        if (p.x > 1) { p.x = 1 - epsilon * rand_uniform();  v.x *= -0.5f; }
        if (p.y < 0) { p.y = epsilon * rand_uniform();  v.y *= -0.5f; }
        if (p.y > 1) { p.y = 1 - epsilon * rand_uniform();  v.y *= -0.5f; }
        if (p.z < 0) { p.z = epsilon * rand_uniform();  v.z *= -0.5f; }
        if (p.z > 1) { p.z = 1 - epsilon * rand_uniform();  v.z *= -0.5f; }

		// update particle
		particles[i].p = p;
	}
}


void simulate_pic(float dt, cgp::numarray<particle_element>& particles, grid_structure& grid, pic_parameters_structure const& pic_parameters) {
    update_grid(particles, grid, pic_parameters, dt);
    update_particle_velocities(particles, grid, pic_parameters, dt);
    update_particle_positions(particles, grid, pic_parameters, dt);
    //std::cout << grid.velocity_y(3, 3, 1) << std::endl;
}