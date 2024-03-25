#include "simulation.hpp"

using namespace cgp;

void solve_incompressibility_cell(grid_structure& grid, int i, int j, int k, pic_parameters_structure pic_parameters, bool rho_correction) {
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
        if (rho_correction) {
            div = div - pic_parameters.stiffness * (grid.rho(i, j, k) - pic_parameters.rho0);
        }
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

void solve_incompressibility(grid_structure& grid, pic_parameters_structure pic_parameters, bool rho_correction) {
    int nb_iter = 10; // gauss-seidel iterations
    for (int iter = 0; iter < nb_iter; ++iter) {
        for (int i = 0; i < grid.cellnumber; ++i) {
            for (int j = 0; j < grid.cellnumber; ++j) {
                for (int k = 0; k < grid.cellnumber; ++k) {
                    solve_incompressibility_cell(grid, i, j, k, pic_parameters, rho_correction);
				}
			}
        }
    }
}

void update_grid(cgp::numarray<particle_element> const& particles, grid_structure& grid, pic_parameters_structure const& pic_parameters, float dt, bool rho_correction) {
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
        float c = grid.cellsize;
        int coords_bottom_left[3] = { floor(part.p[0] / c), floor(part.p[1] / c), floor(part.p[2] / c) };
        float const dx = std::abs((part.p[0] / c) - coords_bottom_left[0] - 0.5f); // -0.5f to be in the middle of the cell
        float const dy = std::abs((part.p[1] / c) - coords_bottom_left[1] - 0.5f);
        float const dz = std::abs((part.p[2] / c) - coords_bottom_left[2] - 0.5f);

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
    solve_incompressibility(grid, pic_parameters, rho_correction);
}

void update_particle_velocities(cgp::numarray<particle_element>& particles, grid_structure const& grid, pic_parameters_structure const& pic_parameters, float dt) {
	// update velocities
	for (int i = 0; i < particles.size(); ++i) {
		// update velocity
        vec3 p = particles[i].p;
        vec3 v = particles[i].v;
		v = vec3(0,0,0);

        int coords_bottom_left[3] = { floor(p[0] / grid.cellsize), floor(p[1] / grid.cellsize), floor(p[2] / grid.cellsize) };
        float bar_coords[3] = { std::abs((p[0] / grid.cellsize) - coords_bottom_left[0] - 0.5f), std::abs((p[1] / grid.cellsize) - coords_bottom_left[1] - 0.5f), std::abs((p[2] / grid.cellsize) - coords_bottom_left[2] - 0.5f)}; // 0.5f to be in the middle of the cell

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
    float h = pic_parameters.h;

    // update positions
	for (int i = 0; i < particles.size(); ++i) {
		vec3 p = particles[i].p;
		vec3 v = particles[i].v;

		// update position
		p += dt * v;

        // Collision
        // small perturbation to avoid alignment
        if (p.x < 0) { p.x = h + epsilon * rand_uniform();  v.x *= -0.5f; }
        if (p.x > 1) { p.x = 1 - h - epsilon * rand_uniform();  v.x *= -0.5f; }
        if (p.y < 0) { p.y = h + epsilon * rand_uniform();  v.y *= -0.5f; }
        if (p.y > 1) { p.y = 1 - h - epsilon * rand_uniform();  v.y *= -0.5f; }
        if (p.z < 0) { p.z = h + epsilon * rand_uniform();  v.z *= -0.5f; }
        if (p.z > 1) { p.z = 1 - h - epsilon * rand_uniform();  v.z *= -0.5f; }

		// update particle
		particles[i].p = p;
	}
}

void update_density_pic(grid_structure& grid, cgp::numarray<particle_element>& particles, float m)
{
    //Compute the density value in each cell of the grid
    grid.rho.fill(0);
    for (int i = 0; i < particles.size(); i++) {
        vec3 p = particles[i].p;
        int coords_bottom_left[3] = { floor(p[0] / grid.cellsize), floor(p[1] / grid.cellsize), floor(p[2] / grid.cellsize) };

        grid.rho(coords_bottom_left[0], coords_bottom_left[1], coords_bottom_left[2]) += m;
    }

    //Normalize the density value by the cell volume
    for (int i = 0; i < grid.cellnumber; i++) {
        for (int j = 0; j < grid.cellnumber; j++) {
            for (int k = 0; k < grid.cellnumber; k++) {
				grid.rho(i, j, k) /= std::pow(grid.cellsize, 3);
			}
		}
	}
}


void simulate_pic(float dt, cgp::numarray<particle_element>& particles, grid_structure& grid, pic_parameters_structure const& pic_parameters, bool rho_correction) {
    update_grid(particles, grid, pic_parameters, dt, rho_correction);
    update_particle_velocities(particles, grid, pic_parameters, dt);
    update_particle_positions(particles, grid, pic_parameters, dt);
    if (rho_correction) {
		update_density_pic(grid, particles, pic_parameters.m);
    }
    update_field(grid, pic_parameters, particles);
}

cgp::grid_3D<cgp::vec3> get_local_average(cgp::numarray<particle_element>& particles, grid_structure const& grid, pic_parameters_structure pic_parameters) {
    float h = pic_parameters.h;

    cgp::int3 samples = grid.implicit_surface.field_param.domain.samples;
    cgp::grid_3D<cgp::vec3> local_avg(samples);
    local_avg.fill(cgp::vec3(0, 0, 0));
    cgp::grid_3D<float> total_value(samples);
    total_value.fill(0);
    cgp::vec3 cell_size = grid.implicit_surface.field_param.domain.voxel_length();
    float R = pic_parameters.R;
    cgp::vec3 nb_cells_to_go = vec3(ceil(R / cell_size.x), ceil(R / cell_size.y), ceil(R / cell_size.z));


    for (int i = 0; i < particles.size(); ++i) {
        vec3 p = particles[i].p;

        // update the local average for the surface reconstruction
        int coords_bottom_left[3] = { floor(p[0] * samples.x), floor(p[1] * samples.y), floor(p[2] * samples.z) };
        for (int ki = coords_bottom_left[0] - nb_cells_to_go.x; ki < coords_bottom_left[0] + nb_cells_to_go.x; ki++) {
            for (int kj = coords_bottom_left[1] - nb_cells_to_go.y; kj < coords_bottom_left[1] + nb_cells_to_go.y; kj++) {
                for (int kk = coords_bottom_left[2] - nb_cells_to_go.z; kk < coords_bottom_left[2] + nb_cells_to_go.z; kk++) {
                    if (ki >= 0 && ki < samples.x && kj >= 0 && kj < samples.y && kk >= 0 && kk < samples.z) {
                        cgp::vec3 position = grid.implicit_surface.field_param.domain.position({ ki, kj, kk });
                        float d = norm(position - p);
                        if (d < R) {
                            float weight = k(d / R);
                            local_avg(ki, kj, kk) += p * weight;
                            total_value(ki, kj, kk) += weight;
                        }
                    }
                }
            }
        }
    }


    // do the average
    for (int ki = 0; ki < samples.x; ki++) {
        for (int kj = 0; kj < samples.y; kj++) {
            for (int kk = 0; kk < samples.z; kk++) {
                if (total_value(ki, kj, kk) > 0) {
                    local_avg(ki, kj, kk) /= total_value(ki, kj, kk);
                }
            }
        }
    }

    return local_avg;
}


void update_field(grid_structure& grid, pic_parameters_structure const& pic_parameters, cgp::numarray<particle_element>& particles)
{
    grid.implicit_surface.field_param.field.fill(0.f);
    float R = pic_parameters.R;
    float r = pic_parameters.radius;

    cgp::spatial_domain_grid_3D domain = grid.implicit_surface.field_param.domain;

    cgp::grid_3D<cgp::vec3> x_bar = get_local_average(particles, grid, pic_parameters);

    // Fill the discrete field values
    for (int kz = 0; kz < domain.samples.z; kz++) {
        for (int ky = 0; ky < domain.samples.y; ky++) {
            for (int kx = 0; kx < domain.samples.x; kx++) {

                vec3 const p = domain.position({ kx, ky, kz });
                grid.implicit_surface.field_param.field(kx, ky, kz) = - (cgp::norm(p - x_bar(kx, ky, kz)) - r); // to have the normal in the correct direction

            }
        }
    }

    // Update the implicit surface
    grid.implicit_surface.update_field_attributes(0.f);
}

float k(float s) {
    return std::max(0.f, (float)pow((1.0f - pow(s, 2)), 3));
}