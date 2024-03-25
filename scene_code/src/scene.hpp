#pragma once


#include "cgp/cgp.hpp"
#include "environment.hpp"

#include "simulation/simulation.hpp"

using cgp::mesh_drawable;


struct gui_parameters {
	bool display_color = false;
	bool display_particles = false;
	bool display_walls = true;
	bool display_grid = false;
	bool rho_correction = true;
	bool display_reconstructed_surface = true;
	float angle_coords = 0.375f;
};

// The structure of the custom scene
struct scene_structure : cgp::scene_inputs_generic {
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //
	/*camera_controller_orbit camera_control;
	camera_projection_orthographic camera_projection;*/
	camera_controller_orbit_euler camera_control;
	camera_projection_perspective camera_projection;
	window_structure window;

	mesh_drawable global_frame;          // The standard global frame
	environment_structure environment;   // Standard environment controler
	input_devices inputs;                // Storage for inputs status (mouse, keyboard, window dimension)
	gui_parameters gui;                  // Standard GUI element storage
	
	// ****************************** //
	// Elements and shapes of the scene
	// ****************************** //
	cgp::timer_basic timer;

	pic_parameters_structure pic_parameters; // Physical parameter related to PIC
	grid_structure grid;
	cgp::numarray<particle_element> particles;      // Storage of the particles
	cgp::mesh_drawable sphere_particle; // Sphere used to display a particle

	cgp::mesh_drawable field_quad; // quad for the ground
	cgp::mesh_drawable wall_quad; // quad for the walls
	cgp::mesh_drawable grid_mesh; // mesh for the grid


	// ****************************** //
	// Functions
	// ****************************** //

	void initialize();    // Standard initialization to be called before the animation loop
	void display_frame(); // The frame display to be called within the animation loop
	void display_gui();   // The display of the GUI, also called within the animation loop

	void initialize_pic();

	void mouse_move_event();
	void mouse_click_event();
	void keyboard_event();
	void idle_frame();

};





