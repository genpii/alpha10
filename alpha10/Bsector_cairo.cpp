#include "stdafx.h"
#include <cairo-win32.h>


using namespace std;

void cairo(const vector<vector<float>>& env, float dangle)
{
	// calculate max value for normalization
	int line = static_cast<int>(distance(env.begin(), env.end()));
	int sample = static_cast<int>(distance(env[0].begin(), env[0].end()));
	vector<float> maxvec(line, 0);
	vector<int>::iterator it;
	int tmp;
	for (int i = 0; i < line; ++i){
		maxvec[i] = *max_element(env[i].begin() + 128, env[i].end() - 128);
	}
	float max = *max_element(maxvec.begin(), maxvec.end());
	
	// draw setting
	double xc = static_cast<double>(sample / 4);
	double yc = 10.0;
	double sangle = (90.0 - static_cast<double>((line * dangle) / 2.0)) * (M_PI / 180.0);
	double eangle;
	double radius;
	double c;
	double gain = 60.0;
	
	cairo_surface_t *surface;
	cairo_t *cr;

	surface = cairo_image_surface_create(CAIRO_FORMAT_RGB24, sample / 2 , sample / 2);
	cr = cairo_create(surface);
	cairo_set_line_width(cr, 1.0);

	//cairo_select_font_face(cr, "serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
	//cairo_set_font_size(cr, 64.0);
	//cairo_set_source_rgb(cr, 0.5, 0.5, 0.5);
	//cairo_move_to(cr, 10.0, 50.0);
	//cairo_show_text(cr, "Gen Onodera");
	for (int i = 0; i < line; ++i){
		eangle = sangle + static_cast<double>(dangle * (M_PI / 180));
		for (int j = 0; j < sample / 4; ++j){
			radius = static_cast<double>(sample / 4 - j);

			c = (20 * log10(env[i][(sample - 1) - 4 * j] / max) + gain) / gain;
			if (c > 1.0) c = 1.0;
			if (c < 0) c = 0;
			cairo_set_source_rgb(cr, c, c, c);
			cairo_arc(cr, xc, yc, radius, sangle, eangle);
			cairo_stroke(cr);
		}
		sangle = eangle;
	}

	cairo_destroy(cr);
	cairo_surface_write_to_png(surface, "hello.png");
	cairo_surface_destroy(surface);


}