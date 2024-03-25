#pragma once

#define _USE_MATH_DEFINES

#include "cgp/cgp.hpp"
#include <math.h> // for M_PI
#include "../display_surface/implicit_surface.hpp"



// PIC Particle
struct particle_element
{
    cgp::vec3 p; // Position
    cgp::vec3 v; // Speed

    particle_element() : p{0,0,0},v{0,0,0} {}
};


enum empty { emptycell, fullcell};

// PIC simulation parameters
struct pic_parameters_structure
{
    // Influence distance of a particle
    float h = 0.008f;

    // Radius used for displaying the particles (larger than the actual radius h so that it is more visible and there is no holes between particles)
    float radius = h * 2.f;

    // R for the surface reconstruction
    float R = radius * 2.f;

    // Rest density (normalized to 1 - real unit should be 1000kg/m^2)
    float rho0 = 1;

    // Total mass of a particle
    float m = rho0 * 4./3. * M_PI * h * h * h / 8.; // we want 8 particles per cell

    // Stiffness to prevent the density to increase too much
    float stiffness = 0.1f;
};

// Grid structure
class grid_structure
{
public:
	cgp::vec3 base[4]; // Square basis of the grid
	float height; // Height of the grid
    int cellnumber; // Number of cells in each direction
	float cellsize; // Size of a cell
    cgp::grid_3D<empty> empty; // status of each cell: empty or not
    cgp::grid_3D<float> velocity_x;
    cgp::grid_3D<float> velocity_y;
    cgp::grid_3D<float> velocity_z;
    cgp::grid_3D<float> rho;
    implicit_surface_structure implicit_surface; // Structures used for the implicit surface

    grid_structure() : base{ cgp::vec3(0,0,0), cgp::vec3(1,0,0), cgp::vec3(1,0,1), cgp::vec3(0,0,1)}, height{1.0f}, cellnumber{30}, cellsize{1.0f / cellnumber} {
        int cellnumber = 30;
        empty.resize(cellnumber, cellnumber, cellnumber);
        empty.fill(emptycell);
        velocity_x.resize(cellnumber, cellnumber, cellnumber);
        velocity_x.fill(0.f);
        velocity_y.resize(cellnumber, cellnumber, cellnumber);
        velocity_y.fill(0.f);
        velocity_z.resize(cellnumber, cellnumber, cellnumber);
        velocity_z.fill(0.f);
        rho.resize(cellnumber, cellnumber, cellnumber);
        rho.fill(0.f);
        // Number of voxels
        int samples = 50;
        // Dimension of the domain
        cgp::vec3 length = { 1,1,1 };
        implicit_surface.set_domain(samples, length, cgp::vec3(0.5,0.5,0.5));
        implicit_surface.field_param.field.resize(samples);
        implicit_surface.field_param.field.fill(0.f);
    }
};

void solve_incompressibility_cell(grid_structure& grid, int i, int j, int k, pic_parameters_structure pic_parameters, bool rho_correction);
void solve_incompressibility(grid_structure& grid, pic_parameters_structure pic_parameters, bool rho_correction);
void update_density_pic(grid_structure& grid, cgp::numarray<particle_element>& particles, float m);
void update_grid(cgp::numarray<particle_element> const& particles, grid_structure& grid, pic_parameters_structure const& sph_parameters, float dt, bool rho_correction);
void update_particle_velocities(cgp::numarray<particle_element>& particles, grid_structure const& grid, pic_parameters_structure const& sph_parameters, float dt);
void update_particle_positions(cgp::numarray<particle_element>& particles, grid_structure& grid, pic_parameters_structure const& sph_parameters, float dt);
void simulate_pic(float dt, cgp::numarray<particle_element>& particles, grid_structure& grid, pic_parameters_structure const& sph_parameters, bool rho_correction);
cgp::grid_3D<cgp::vec3> get_local_average(cgp::numarray<particle_element>& particles, grid_structure const& grid, pic_parameters_structure pic_parameters);
void update_field(grid_structure& grid, pic_parameters_structure const& pic_parameters, cgp::numarray<particle_element>& particles);
float k(float s);
