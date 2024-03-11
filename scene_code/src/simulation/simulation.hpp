#pragma once

#include "cgp/cgp.hpp"



// SPH Particle
struct particle_element
{
    cgp::vec3 p; // Position
    cgp::vec3 v; // Speed
    cgp::vec3 f; // Force

    float rho;      // density at this particle position
    float pressure; // pressure at this particle position

    particle_element() : p{0,0,0},v{0,0,0},f{0,0,0},rho(0),pressure(0) {}
};

// SPH simulation parameters
struct sph_parameters_structure
{
    // Influence distance of a particle (size of the kernel)
    float h = 0.12f;

    // Rest density (normalized to 1 - real unit should be 1000kg/m^2)
    float rho0 = 1;

     // Total mass of a particle (consider rho0 h^2)
    float m = rho0*h*h;

    // viscosity parameter
    float nu = 0.02f;   
     
    // Stiffness converting density to pressure
    float stiffness = 8.0f;
    
};


void simulate(float dt, cgp::numarray<particle_element>& particles, sph_parameters_structure const& sph_parameters);

enum status { rigid, solid};
enum empty { emptycell, fullcell};

// PIC simulation parameters
struct pic_parameters_structure
{
    // Influence distance of a particle (size of the kernel)
    float h = 0.03f;

    // Rest density (normalized to 1 - real unit should be 1000kg/m^2)
    float rho0 = 1;

    // Total mass of a particle (consider rho0 h^2)
    float m = rho0 * h * h;

};

// Grid structure
class grid_structure
{
public:
	cgp::vec3 base[4]; // Square basis of the grid
	float height; // Height of the grid
    int cellnumber; // Number of cells in each direction
	float cellsize; // Size of a cell
    cgp::grid_3D<status> status; // status of each cell: rigid movement or fluid
    cgp::grid_3D<empty> empty; // status of each cell: empty or not
    cgp::grid_3D<cgp::mat3> stress; // stress matrix for each cell
    cgp::grid_3D<float> velocity_x;
    cgp::grid_3D<float> velocity_y;
    cgp::grid_3D<float> velocity_z;

    grid_structure() : base{ cgp::vec3(0,0,0), cgp::vec3(1,0,0), cgp::vec3(1,0,1), cgp::vec3(0,0,1)}, height{1.0f}, cellnumber{15}, cellsize{1.0f / cellnumber} {
        int cellnumber = 15;
        status.resize(cellnumber, cellnumber, cellnumber);
        status.fill(rigid);
        empty.resize(cellnumber, cellnumber, cellnumber);
        empty.fill(emptycell);
        stress.resize(cellnumber, cellnumber, cellnumber);
        stress.fill(cgp::mat3(0.f));
        velocity_x.resize(cellnumber, cellnumber, cellnumber);
        velocity_x.fill(0.f);
        velocity_y.resize(cellnumber, cellnumber, cellnumber);
        velocity_y.fill(0.f);
        velocity_z.resize(cellnumber, cellnumber, cellnumber);
        velocity_z.fill(0.f);
    }
};

void solve_incompressibility_cell(grid_structure& grid, int i, int j, int k);
void solve_incompressibility(grid_structure& grid);
void update_grid(cgp::numarray<particle_element> const& particles, grid_structure& grid, pic_parameters_structure const& sph_parameters, float dt);
void update_particle_velocities(cgp::numarray<particle_element>& particles, grid_structure const& grid, pic_parameters_structure const& sph_parameters, float dt);
void update_particle_positions(cgp::numarray<particle_element>& particles, grid_structure& grid, pic_parameters_structure const& sph_parameters, float dt);
void simulate_pic(float dt, cgp::numarray<particle_element>& particles, grid_structure& grid, pic_parameters_structure const& sph_parameters);