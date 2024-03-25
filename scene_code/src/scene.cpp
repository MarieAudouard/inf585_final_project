#include "scene.hpp"
#include "display_surface/implicit_surface.hpp"


using namespace cgp;

void scene_structure::initialize()
{
	std::cout << "Initialization start start" << std::endl;
	camera_control.initialize(inputs, window); // Give access to the inputs and window global state to the camera controler
	camera_control.look_at({ 0.0f, 0.0f, 2.0f }, {0,0,0}, {0,1,0});
	global_frame.initialize_data_on_gpu(mesh_primitive_frame());

	// bottom points of the grid
	vec3 p000 = vec3(0, 0, 0);
	vec3 p100 = vec3(1, 0, 0);
	vec3 p110 = vec3(1, 0, 1);
	vec3 p010 = vec3(0, 0, 1);

	// set the grid
	grid.base[0] = p000;
	grid.base[1] = p100;
	grid.base[2] = p110;
	grid.base[3] = p010;

	//height of the grid
	float const h = 1.0f;

	// set the height of the grid
	grid.height = h;

	// for the ground
	field_quad.initialize_data_on_gpu(mesh_primitive_quadrangle(p010, p110, p100, p000) );
	field_quad.material.color = { 0.95f, 0.95f, 0.95f };
	field_quad.material.phong.specular = 0.0f; // non-specular surface

	// for the walls
	wall_quad.initialize_data_on_gpu(mesh_primitive_quadrangle(p000 + vec3(0,h,0), p010 + vec3(0, h, 0), p010, p000));
	wall_quad.material.color = { 0.99f, 0.99f, 0.99f };
	wall_quad.material.phong.specular = 0.0f; // non-specular surface

	std::cout << "Initialization start" << std::endl;
	initialize_pic();
	sphere_particle.initialize_data_on_gpu(mesh_primitive_sphere(1.0,{0,0,0},10,10));
	sphere_particle.model.scaling = pic_parameters.radius;
	sphere_particle.material.color = { 1,1,0 }; // yellow

	// for the grid
	grid_mesh.initialize_data_on_gpu(mesh_primitive_cubic_grid(p010 + vec3(0, grid.cellsize/2, 0), p110 + vec3(0, grid.cellsize / 2, 0), p100 + vec3(0, grid.cellsize / 2, 0), p000 + vec3(0, grid.cellsize / 2, 0), p010 + vec3(0, h + grid.cellsize / 2, 0), p110 + vec3(0, h + grid.cellsize / 2, 0), p100 + vec3(0, h + grid.cellsize / 2, 0), p000 + vec3(0, h + grid.cellsize / 2, 0), grid.cellnumber, grid.cellnumber, grid.cellnumber));
	grid_mesh.material.color = { 0.7f, 0.f, 0.f };
	std::cout << "Initialization done" << std::endl;
}

void scene_structure::initialize_pic() {
	// Initial particle spacing (8 particles per cell according to the paper)
	float const c = grid.cellsize / 2;

	// Fill the column with an angle in (angle_coords,0,angle_coords) and with 1/4 width
	particles.clear();

	for (float x = gui.angle_coords + c / 2; x < gui.angle_coords + 0.25f; x = x + c)
	{
		for (float y = c / 2; y < 0.5f; y = y + c)
		{
			for (float z = gui.angle_coords + c / 2; z < gui.angle_coords + 0.25f; z = z + c)
			{
				particle_element particle;
				particle.p = { x + c / 8.0 * rand_uniform(),y + c / 8.0 * rand_uniform(),z + c / 8.0 * rand_uniform() };
				particles.push_back(particle);
				int coords_bottom_left[3] = { floor(particle.p[0] / grid.cellsize), floor(particle.p[1] / grid.cellsize), floor(particle.p[2] / grid.cellsize) };
				grid.empty(coords_bottom_left[0], coords_bottom_left[1], coords_bottom_left[2]) = fullcell;
			}
		}
	}

	update_density_pic(grid, particles, pic_parameters.m);
	update_field(grid, pic_parameters, particles);
	grid.implicit_surface.update_field_attributes(0.f);
}

void scene_structure::display_frame()
{
	// Set the light to the current position of the camera
	environment.light = camera_control.camera_model.position();
	
	timer.update(); // update the timer to the current elapsed time
	float const dt = 0.003f * timer.scale;
	simulate_pic(dt, particles, grid, pic_parameters, gui.rho_correction);

	if (gui.display_walls) {
		draw(field_quad, environment);
		//draw the walls that are opposite of the camera using transformations
		wall_quad.model.translation = { 0, 1, 0 };
		wall_quad.model.rotation = rotation_transform::from_axis_angle({ 1, 0, 0 }, 3.1415f / 2);
		if (dot(camera_control.camera_model.center_of_rotation - camera_control.camera_model.position(), vec3(-1,0,0)) > 0)
			draw(wall_quad, environment);
		wall_quad.model.translation = { 1, 0, 0 };
		wall_quad.model.rotation = rotation_transform::from_axis_angle({ 0, 1, 0 }, - 3.1415f / 2);
		if (dot(camera_control.camera_model.center_of_rotation - camera_control.camera_model.position(), vec3(0, 0, -1)) > 0)
			draw(wall_quad, environment);
		wall_quad.model.translation = { 1, 0, 1};
		wall_quad.model.rotation = rotation_transform::from_axis_angle({ 0, 1, 0 }, 3.1415f);
		if (dot(camera_control.camera_model.center_of_rotation - camera_control.camera_model.position(), vec3(1, 0, 0)) > 0)
			draw(wall_quad, environment);
		wall_quad.model.translation = { 0, 0, 1 };
		wall_quad.model.rotation = rotation_transform::from_axis_angle({ 0, 1, 0 }, 3.1415f / 2);
		if (dot(camera_control.camera_model.center_of_rotation - camera_control.camera_model.position(), vec3(0, 0, 1)) > 0)
			draw(wall_quad, environment);
	}

	if (gui.display_particles && !gui.display_color) {
		for (int k = 0; k < particles.size(); ++k) {
			vec3 const& p = particles[k].p;
			sphere_particle.model.translation = p;
			sphere_particle.model.scaling = pic_parameters.h;
			draw(sphere_particle, environment);
		}
	}

	if (gui.display_color) {
		for (int k = 0; k < particles.size(); ++k) {
			vec3 const& p = particles[k].p;
			sphere_particle.model.translation = p;
			sphere_particle.model.scaling = pic_parameters.h + 0.008f;
			draw(sphere_particle, environment);
		}
	}

	if (gui.display_grid) {
		draw_wireframe(grid_mesh, environment);
	}

	if (gui.display_reconstructed_surface) {
		grid.implicit_surface.update_field_attributes(0.f);
		draw(grid.implicit_surface.drawable_param.shape, environment);
	}

}

void scene_structure::display_gui()
{
	ImGui::SliderFloat("Timer scale", &timer.scale, 0.01f, 4.0f, "%0.2f");
	ImGui::SliderFloat("Position of the angle", &gui.angle_coords, 0.0f, 0.75f);

	bool const restart = ImGui::Button("Restart");
	if (restart)
		initialize_pic();

	ImGui::Checkbox("Particles for display", &gui.display_color);
	ImGui::Checkbox("Real Particles", &gui.display_particles);
	ImGui::Checkbox("Display walls", &gui.display_walls);
	ImGui::Checkbox("Display grid", &gui.display_grid);
	ImGui::Checkbox("Display reconstructed surface", &gui.display_reconstructed_surface);
	ImGui::Checkbox("Correct particle position using density", &gui.rho_correction);
}

void scene_structure::mouse_move_event()
{
	if (!inputs.keyboard.shift)
		camera_control.action_mouse_move(environment.camera_view);
}
void scene_structure::mouse_click_event()
{
	camera_control.action_mouse_click(environment.camera_view);
}
void scene_structure::keyboard_event()
{
	camera_control.action_keyboard(environment.camera_view);
}
void scene_structure::idle_frame()
{
	camera_control.idle_frame(environment.camera_view);
}

