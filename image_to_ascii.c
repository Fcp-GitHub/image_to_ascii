/*
 *	image_to_ascii.c -- convert a given image to an ASCII matrix and print it.	
 *
 *	Copyright (C) 2024  fcp.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  ADDITIONAL TERMS: Preservation of author attributions is required. 
 *  Any changes to the program must be documented, as stated in the GNU GPLv3 
 *  itself. A change notice should contain the name of the author of said
 *  modifications, the date when they are made public and a short explanation
 *  of what has been done.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>
#ifdef __GNUC__
#include <getopt.h>	// For getopt_long().
#endif
#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>	// For image loading.
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include <stb_image_resize2.h>	// For image resizing.

#define TERM_WIDTH 80 
#define TERM_HEIGHT 24
#define ERROR_CHAR 0
#define ASCII_PRECISION 11	// Number of ASCII characters used to convert the image.
#define COLOR_PRECISION 8	// Number of different colors that are supported by this implementation.
#define EPSILON 0.00001		// Kind of "color epsilon" used to detect black and white pixels.

// Quick macro to ensure that the resulting hue value will be correct (used in get_hsl()).
#define POLAR_HUE(N) (N) >= 0 ? (N) : (360 + (N))
					
// Enumeration used to keep track of different errors.
typedef enum {
	ALL_GOOD,
	IMAGE_LOADING_ERROR,
	IMAGE_RESIZING_ERROR,
	CHOOSE_ASCII_ERROR,
	USAGE_ERROR
} error_t;

// ASCII characters lookup table.
char ascii_table[ASCII_PRECISION] = {
	' ',
	'.',
	',',
	':',
	';',
	't',
	'x',	
	'@',
	'#',	
	'N',
	'M'	
};

// struct used to manage rgb colors more easily.
typedef struct {
	double r, g, b;
} rgb_color;

// 3-bit RGB colors lookup table.
rgb_color color_table[COLOR_PRECISION] = {
	{0,0,0},	// black 	(ANSI code 30)
	{1,0,0},	// red	 	(ANSI code 31)
	{0,1,0},	// green 	(ANSI code 32)
	{1,1,0},	// yellow	(ANSI code 33)
	{0,0,1},	// blue		(ANSI code 34)
	{1,0,1},	// magenta	(ANSI code 35)
	{0,1,1},	// cyan		(ANSI code 36)
	{1,1,1}		// white	(ANSI code 37)
};

// struct used to manage hsl colors more easily.
typedef struct {
	double h, s, l;
} hsl_color;

// HSL values of the above 3-bit colors.
hsl_color hsl_table[COLOR_PRECISION] = {
	{0,0,0},	// black
	{0,1,.5},	// red
	{120,1,.5},	// green
	{60,1,.5},	// yellow
	{240,1,.5},	// blue
	{300,1,.1},	// magenta
	{180,1,.1},	// cyan
	{0,0,1}		// white
};

// Output functions.
//-------------------------------------------------------------------------------
// Print error explanations.
void print_error(error_t error);
// Print help message.
void print_help(void);
// Print char matrix.
void print_ascii(const char result[TERM_WIDTH][TERM_HEIGHT]);
// Print colored char matrix.
void print_colored_ascii(const char result[TERM_WIDTH][TERM_HEIGHT], int color_map[TERM_WIDTH][TERM_HEIGHT], const char choice);

// Image loading and resizing
//-------------------------------------------------------------------------------
error_t get_image_data(const char* filename, unsigned char** data, int* nrChannels);

// Required math calculations.
//-------------------------------------------------------------------------------
// Luminance.
double Y(const rgb_color c);
// Perceived lightness.
double L_star(const double Y);
// Get hue from rgb representation.
const hsl_color get_hsl(const rgb_color c);

// ASCII matrix generation.
//-------------------------------------------------------------------------------
// Choose ASCII character based on perceived lightness.
const char choose_ascii(double L_star);
// Find the most suitable color from the hue of the analyzed pixel.
const int choose_color(const rgb_color c);
// Generate black and white ASCII matrix based on perceived lightness.
error_t grayscale_ascii(unsigned char* data, char result[TERM_WIDTH][TERM_HEIGHT], int nrChannels);
// Generate colored ASCII matrix based on perceived lightness and color similarity.
error_t colored_ascii(unsigned char* data, char result[TERM_WIDTH][TERM_HEIGHT], int color_map[TERM_WIDTH][TERM_HEIGHT], int nrChannels);

int main(int argc, char* argv[])
{
#ifdef __GNUC__
	// Parse command-line argument.
	// Long options.
	static struct option long_options[] = 
	{
		{"grayscale", required_argument, NULL, 'g'},
		{"colored", required_argument, NULL, 'c'},
		{"full", required_argument, NULL, 'f'},
		{"negative", required_argument, NULL, 'n'},
		{"help", no_argument, NULL, 'h'},
		{"usage", no_argument, NULL, 'h'}
	};
	int opt = getopt_long(argc, argv, ":g:c:f:n:hu", long_options, NULL);
#else
	// Short options.
	int opt = getopt(argc, argv, ":g:c:f:n:hu");	// Store the argument.
#endif
	char ascii_map[TERM_WIDTH][TERM_HEIGHT] = {0};	// the ASCII map will be stored here. 
	int color_map[TERM_WIDTH][TERM_HEIGHT] = {0};	// color_map that will be filled if colored methods are required.
	unsigned char* data = {NULL};	// file data.
	int nrChannels = 0;	// number of color channels used.
	error_t error;	// Any potential error will be stored here.

	switch(opt)
	{
		case 'g' :	// Grayscale ASCII.
			if ((error = get_image_data(optarg, &data, &nrChannels)) != ALL_GOOD)
			{
				print_error(error);
				return -1;
			}
			if ((error = grayscale_ascii(data, ascii_map, nrChannels)) != ALL_GOOD)
			{
				print_error(error);
				return -1;
			}
			// Print final result.
			print_ascii(ascii_map);
			break;
		case 'c' :	// Colored ASCII. 
			if ((error = get_image_data(optarg, &data, &nrChannels)) != ALL_GOOD)
			{
				print_error(error);
				return -1;
			}
			if ((error = colored_ascii(data, ascii_map, color_map, nrChannels)) != ALL_GOOD)
			{
				print_error(error);
				return -1;
			}
			// Print final result.
			print_colored_ascii(ascii_map, color_map, 'c');
			break;
		case 'f' :	// Full colored ASCII. 
			if ((error = get_image_data(optarg, &data, &nrChannels)) != ALL_GOOD)
			{
				print_error(error);
				return -1;
			}
			if ((error = colored_ascii(data, ascii_map, color_map, nrChannels)) != ALL_GOOD)
			{
				print_error(error);
				return -1;
			}
			// Print final result.
			print_colored_ascii(ascii_map, color_map, 'f');
			break;
		case 'n' :	// Full colored "negative" ASCII. 
			if ((error = get_image_data(optarg, &data, &nrChannels)) != ALL_GOOD)
			{
				print_error(error);
				return -1;
			}
			if ((error = colored_ascii(data, ascii_map, color_map, nrChannels)) != ALL_GOOD)
			{
				print_error(error);
				return -1;
			}
			// Print final result.
			print_colored_ascii(ascii_map, color_map, 'n');
			break;
		case 'h' : print_help(); return 0;
		case 'u' : print_error(USAGE_ERROR); return 0;
		case ':' : print_error(USAGE_ERROR); return -1;
		default  : print_error(USAGE_ERROR); return -1;
	}
		
	// Free image's data.
	stbi_image_free(data);

	return 0;
}

// Output functions.
//-------------------------------------------------------------------------------
void print_error(error_t error)
{
	switch(error)
	{
		case IMAGE_LOADING_ERROR:
			puts("Error during image loading.");
			break;
		case IMAGE_RESIZING_ERROR:
			puts("Error during image resizing.");
			break;
		case CHOOSE_ASCII_ERROR:
			puts("Error finding the right ASCII character.");
			break;
		case USAGE_ERROR:
			puts("image_to_ascii Copyright (C) 2024 fcp");
			puts("Usage: image_to_ascii [OPTION] [IMAGE]");
			break;
		default:
			break;
	}
}

void print_help(void)
{
	puts("Usage: image_to_ascii [OPTION] [IMAGE]");
	printf("\t-g, --grayscale %10clack and white ASCII characters\n", 'b');
	printf("\t-c, --colored %12colored ASCII characters\n", 'c');
	printf("\t-f, --full-color %9cSCII characters with same background and foreground color\n", 'A');
	printf("\t-n, --negative %11clack ASCII characters with background color\n", 'b');
	printf("\t-h, --help %15cisplay this list\n", 't');
	printf("\t-u, --usage %14cive a short usage message\n", 'g');
	puts("");
}

void print_ascii(const char result[TERM_WIDTH][TERM_HEIGHT])
{
	for (size_t i = 0; i < TERM_HEIGHT; i++)
	{
		for (size_t j = 0; j < TERM_WIDTH; j++)
		{
			printf("%c", result[j][i]);
		}
		puts("");
	}		
}

void print_colored_ascii(const char ascii_map[TERM_WIDTH][TERM_HEIGHT], int color_map[TERM_WIDTH][TERM_HEIGHT], const char choice)
{
	switch(choice)
	{
		case 'c': 
			for (size_t i = 0; i < TERM_HEIGHT; i++)
			{
				for (size_t j = 0; j < TERM_WIDTH; j++)
				{
					printf("\x1b[%dm%c%s", color_map[j][i], ascii_map[j][i], "\x1b[m");
				}
				puts("");
			}
	   		break;
   		case 'f':
			for (size_t i = 0; i < TERM_HEIGHT; i++)
			{	
				for (size_t j = 0; j < TERM_WIDTH; j++)
				{
					printf("\x1b[%d;%dm%c%s", color_map[j][i], color_map[j][i]+10, ascii_map[j][i], "\x1b[m");
				}
				puts("");
			}
			break;
		case 'n':
			for (size_t i = 0; i < TERM_HEIGHT; i++)
			{
				for (size_t j = 0; j < TERM_WIDTH; j++)
				{
					printf("\x1b[%d;%dm%c%s", 30, color_map[j][i]+10, ascii_map[j][i], "\x1b[m");
				}
				puts("");
			}
			break;
	}
}

// Image loading and resizing
//-------------------------------------------------------------------------------
error_t get_image_data(const char* filename, unsigned char** data, int* nrChannels)
{
	// Load the image and resize it.
	int width, height;

	*data = stbi_load(filename, &width, &height, nrChannels, 0);
	
	// handle loading error.
	if (!(*data)) { stbi_image_free(*data); return IMAGE_LOADING_ERROR; }

	*data = stbir_resize_uint8_linear(*data, width, height, 0, NULL, TERM_WIDTH, TERM_HEIGHT, 0, (stbir_pixel_layout)(*nrChannels));
	// handle resizing error.
	if (!(*data)) { stbi_image_free(*data); return IMAGE_RESIZING_ERROR; }

	return ALL_GOOD;
}

// Required math calculations.
//-------------------------------------------------------------------------------
double Y(const rgb_color c)
{
	return c.r * 0.2126 + c.g * 0.7152 + c.b * 0.0722;
}

double L_star(const double Y)
{
	if (Y <= 0.008856)
		return Y * 903.3;
	else
		return cbrt(Y) * 116 - 16;	
}

// Comparison function used by qsort.
int comparison_function(const void* e1, const void* e2)
{
	double first  = *(const double*)e1;
	double second = *(const double*)e2;
	return (first > second) - (first < second);	
}

const hsl_color get_hsl(const rgb_color c)
{
	// From https://en.wikipedia.org/wiki/HSL_and_HSV#From_RGB.
	// 1. Get max and min rgb components (sort a double[3] array in order to quickly get them).
	double temp[3] = {c.r, c.g, c.b};
	qsort(temp, 3, sizeof(double), comparison_function);
	double x_min = temp[0], x_max = temp[2]; 

	// 2. Get chroma value (for convenience) and actual "value" value.
	double chroma = x_max - x_min;
	double lightness = (x_max + x_min)/2;

	// 3. Get right hue value.
	double hue = 0.f;
	if (chroma != 0.)
	{
		if (x_max == c.r)
			hue = POLAR_HUE(60 * fmod(( (c.g - c.b) / chroma), 6));
		else if (x_max == c.g)
			hue = POLAR_HUE(60 * ( ( (c.b - c.r) / chroma) + 2));
		else 
			hue = POLAR_HUE(60 * ( ( (c.r - c.g) / chroma) + 4));
	}

	// 4. Get saturation value.
	double saturation = 0.f;
	if (lightness != 0.f && lightness != 1.f)
		saturation = fabs( chroma / (1 - fabs(2*x_max - chroma - 1) ) );

	return (hsl_color){hue, saturation, lightness};
}

// ASCII matrix generation.
//-------------------------------------------------------------------------------
const char choose_ascii(const double L_star)
{
	// Iterate through all the lightness levels with a generic for loop,
	// allowing to easily modify the ASCII precision.
	double subinterval = 100.f / (double)(ASCII_PRECISION - 1);
	for (int i = 0; i < ASCII_PRECISION; i++)
	{
		if (L_star <= (i * subinterval))
		{
			return ascii_table[i];
		}	
	}

	// Signal if something impossible (the algorithm shouldn't get here) has happened. 
	return ERROR_CHAR;
}

const int choose_color(const rgb_color c)
{
	// Save current smallest color difference and the required ANSI code. 
	double temp = INFINITY;	 	
	int ansi_code = 29;	
	const hsl_color c1 = get_hsl(c);
	if (c1.l <= EPSILON)	// black
		return ansi_code + 1;
	if (c1.l >= 1-EPSILON)	// white
		return ansi_code + COLOR_PRECISION;
	
	// Iterate through all the available colors and find the best match.
	ansi_code++;
	for (size_t i = 1; i < COLOR_PRECISION-1; i++)
	{
		const double c_diff = fabs( c1.h - hsl_table[i].h );
		if (c_diff < temp)
		{
			temp = c_diff;
			ansi_code++;
		}
	}

	return ansi_code;
}

error_t grayscale_ascii(unsigned char* data, char ascii_map[TERM_WIDTH][TERM_HEIGHT], int nrChannels)
{
	// Iterate through the image's pixels.
	unsigned int byte_per_pixel = nrChannels;
	for(size_t i = 0; i < TERM_WIDTH; i++)
	{
		for(size_t j = 0; j < TERM_HEIGHT; j++)
		{
			// Find correct pixel.
			unsigned char* pixel_offset = data + (i + TERM_WIDTH * j) * byte_per_pixel;
			double r = pixel_offset[0], g = pixel_offset[1], b = pixel_offset[2],
		   				  a = nrChannels >= 4 ? pixel_offset[3] : 0xff;

			// Normalize pixel values.
			r /= 255.; g /= 255.; b /= 255.;

			// Compute perceived lightness.
			double Ls = L_star( Y((rgb_color){r,g,b}) );

			// Find correct character and assign it to result.
			ascii_map[i][j] = choose_ascii(Ls);
			if (ascii_map[i][j] == ERROR_CHAR)	return CHOOSE_ASCII_ERROR;
		}
	}

	return ALL_GOOD;
}

error_t colored_ascii(unsigned char* data, char ascii_map[TERM_WIDTH][TERM_HEIGHT], int color_map[TERM_WIDTH][TERM_HEIGHT], int nrChannels)
{
	// Iterate through the image's pixels.
	unsigned int byte_per_pixel = nrChannels;
	for(size_t i = 0; i < TERM_WIDTH; i++)
	{
		for(size_t j = 0; j < TERM_HEIGHT; j++)
		{
			// Find correct pixel.
			unsigned char* pixel_offset = data + (i + TERM_WIDTH * j) * byte_per_pixel;
			double r = pixel_offset[0], g = pixel_offset[1], b = pixel_offset[2],
		   				  a = nrChannels >= 4 ? pixel_offset[3] : 0xff;

			// Normalize pixel values.
			r /= 255.; g /= 255.; b /= 255.;

			// Compute perceived lightness.
			double Ls = L_star( Y((rgb_color){r,g,b}) );

			// Find correct character and assign it to ascii_map.
			ascii_map[i][j] = choose_ascii(Ls);
			if (ascii_map[i][j] == ERROR_CHAR)	return CHOOSE_ASCII_ERROR;

			// Find right color.
			color_map[i][j] = choose_color((rgb_color){r,g,b});
		}
	}

	return ALL_GOOD;
}
