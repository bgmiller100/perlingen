
/*
 * Octave2D.cpp
 *
 * Creates an 'octave noise' pattern (multiple octaves of basic Perlin
 * noise, evaluated at a set of points provided by the user. A grid with
 * unit spacing is used for the Perlin grid - input points should be
 * transformed appropriately to generate features of desired size and 
 * orientation. The user may also specify the number of octaves, and the
 * factor successively applied to each octave. This seeded version of the
 * code also takes as input a permutation table of the numbers 0 to 255
 * and an m x 2 set of offsets in 2D space to apply to the Perlin grid for
 * each octave. m must therefore be at least as large as N_freq.
 *
 * Calling syntax in MATLAB is
 *
 *		noisefield = Octave2D(points, N_freq, fade_factor, Ps, offset)
 *
 * This is a MEX file for MATLAB.
*/

#include <Python.h>
#include <stdio.h>
#include <math.h>

//printf("methods launch");       
//Define Perlin vectors (256 points distributed evenly around circle)
static double vecs_x[256] = 
    {1.0000, 0.9997, 0.9988, 0.9973, 0.9952, 0.9925, 0.9892, 0.9853, 0.9808,
     0.9757, 0.9700, 0.9638, 0.9569, 0.9495, 0.9415, 0.9330, 0.9239, 0.9142,
     0.9040, 0.8932, 0.8819, 0.8701, 0.8577, 0.8449, 0.8315, 0.8176, 0.8032,
     0.7883, 0.7730, 0.7572, 0.7410, 0.7242, 0.7071, 0.6895, 0.6716, 0.6532,
     0.6344, 0.6152, 0.5957, 0.5758, 0.5556, 0.5350, 0.5141, 0.4929, 0.4714,
     0.4496, 0.4276, 0.4052, 0.3827, 0.3599, 0.3369, 0.3137, 0.2903, 0.2667,
     0.2430, 0.2191, 0.1951, 0.1710, 0.1467, 0.1224, 0.0980, 0.0736, 0.0491,
     0.0245, 0.0000, -0.0245, -0.0491, -0.0736, -0.0980, -0.1224, -0.1467, -0.1710,
     -0.1951, -0.2191, -0.2430, -0.2667, -0.2903, -0.3137, -0.3369, -0.3599, -0.3827,
     -0.4052, -0.4276, -0.4496, -0.4714, -0.4929, -0.5141, -0.5350, -0.5556, -0.5758,
     -0.5957, -0.6152, -0.6344, -0.6532, -0.6716, -0.6895, -0.7071, -0.7242, -0.7410,
     -0.7572, -0.7730, -0.7883, -0.8032, -0.8176, -0.8315, -0.8449, -0.8577, -0.8701,
     -0.8819, -0.8932, -0.9040, -0.9142, -0.9239, -0.9330, -0.9415, -0.9495, -0.9569,
     -0.9638, -0.9700, -0.9757, -0.9808, -0.9853, -0.9892, -0.9925, -0.9952, -0.9973,
     -0.9988, -0.9997, -1.0000, -0.9997, -0.9988, -0.9973, -0.9952, -0.9925, -0.9892,
     -0.9853, -0.9808, -0.9757, -0.9700, -0.9638, -0.9569, -0.9495, -0.9415, -0.9330,
     -0.9239, -0.9142, -0.9040, -0.8932, -0.8819, -0.8701, -0.8577, -0.8449, -0.8315,
     -0.8176, -0.8032, -0.7883, -0.7730, -0.7572, -0.7410, -0.7242, -0.7071, -0.6895,
     -0.6716, -0.6532, -0.6344, -0.6152, -0.5957, -0.5758, -0.5556, -0.5350, -0.5141,
     -0.4929, -0.4714, -0.4496, -0.4276, -0.4052, -0.3827, -0.3599, -0.3369, -0.3137,
     -0.2903, -0.2667, -0.2430, -0.2191, -0.1951, -0.1710, -0.1467, -0.1224, -0.0980,
     -0.0736, -0.0491, -0.0245, -0.0000, 0.0245, 0.0491, 0.0736, 0.0980, 0.1224,
     0.1467, 0.1710, 0.1951, 0.2191, 0.2430, 0.2667, 0.2903, 0.3137, 0.3369,
     0.3599, 0.3827, 0.4052, 0.4276, 0.4496, 0.4714, 0.4929, 0.5141, 0.5350,
     0.5556, 0.5758, 0.5957, 0.6152, 0.6344, 0.6532, 0.6716, 0.6895, 0.7071,
     0.7242, 0.7410, 0.7572, 0.7730, 0.7883, 0.8032, 0.8176, 0.8315, 0.8449,
     0.8577, 0.8701, 0.8819, 0.8932, 0.9040, 0.9142, 0.9239, 0.9330, 0.9415,
     0.9495, 0.9569, 0.9638, 0.9700, 0.9757, 0.9808, 0.9853, 0.9892, 0.9925,
     0.9952, 0.9973, 0.9988, 0.9997};
static double vecs_y[256] = 
    {0.0000, 0.0245, 0.0491, 0.0736, 0.0980, 0.1224, 0.1467, 0.1710, 0.1951,
     0.2191, 0.2430, 0.2667, 0.2903, 0.3137, 0.3369, 0.3599, 0.3827, 0.4052,
     0.4276, 0.4496, 0.4714, 0.4929, 0.5141, 0.5350, 0.5556, 0.5758, 0.5957,
     0.6152, 0.6344, 0.6532, 0.6716, 0.6895, 0.7071, 0.7242, 0.7410, 0.7572,
     0.7730, 0.7883, 0.8032, 0.8176, 0.8315, 0.8449, 0.8577, 0.8701, 0.8819,
     0.8932, 0.9040, 0.9142, 0.9239, 0.9330, 0.9415, 0.9495, 0.9569, 0.9638,
     0.9700, 0.9757, 0.9808, 0.9853, 0.9892, 0.9925, 0.9952, 0.9973, 0.9988,
     0.9997, 1.0000, 0.9997, 0.9988, 0.9973, 0.9952, 0.9925, 0.9892, 0.9853,
     0.9808, 0.9757, 0.9700, 0.9638, 0.9569, 0.9495, 0.9415, 0.9330, 0.9239,
     0.9142, 0.9040, 0.8932, 0.8819, 0.8701, 0.8577, 0.8449, 0.8315, 0.8176,
     0.8032, 0.7883, 0.7730, 0.7572, 0.7410, 0.7242, 0.7071, 0.6895, 0.6716,
     0.6532, 0.6344, 0.6152, 0.5957, 0.5758, 0.5556, 0.5350, 0.5141, 0.4929,
     0.4714, 0.4496, 0.4276, 0.4052, 0.3827, 0.3599, 0.3369, 0.3137, 0.2903,
     0.2667, 0.2430, 0.2191, 0.1951, 0.1710, 0.1467, 0.1224, 0.0980, 0.0736,
     0.0491, 0.0245, 0.0000, -0.0245, -0.0491, -0.0736, -0.0980, -0.1224, -0.1467,
     -0.1710, -0.1951, -0.2191, -0.2430, -0.2667, -0.2903, -0.3137, -0.3369, -0.3599,
     -0.3827, -0.4052, -0.4276, -0.4496, -0.4714, -0.4929, -0.5141, -0.5350, -0.5556,
     -0.5758, -0.5957, -0.6152, -0.6344, -0.6532, -0.6716, -0.6895, -0.7071, -0.7242,
     -0.7410, -0.7572, -0.7730, -0.7883, -0.8032, -0.8176, -0.8315, -0.8449, -0.8577,
     -0.8701, -0.8819, -0.8932, -0.9040, -0.9142, -0.9239, -0.9330, -0.9415, -0.9495,
     -0.9569, -0.9638, -0.9700, -0.9757, -0.9808, -0.9853, -0.9892, -0.9925, -0.9952,
     -0.9973, -0.9988, -0.9997, -1.0000, -0.9997, -0.9988, -0.9973, -0.9952, -0.9925,
     -0.9892, -0.9853, -0.9808, -0.9757, -0.9700, -0.9638, -0.9569, -0.9495, -0.9415,
     -0.9330, -0.9239, -0.9142, -0.9040, -0.8932, -0.8819, -0.8701, -0.8577, -0.8449,
     -0.8315, -0.8176, -0.8032, -0.7883, -0.7730, -0.7572, -0.7410, -0.7242, -0.7071,
     -0.6895, -0.6716, -0.6532, -0.6344, -0.6152, -0.5957, -0.5758, -0.5556, -0.5350,
     -0.5141, -0.4929, -0.4714, -0.4496, -0.4276, -0.4052, -0.3827, -0.3599, -0.3369,
     -0.3137, -0.2903, -0.2667, -0.2430, -0.2191, -0.1951, -0.1710, -0.1467, -0.1224,
     -0.0980, -0.0736, -0.0491, -0.0245, 
    };
    
typedef struct Vector2_s {
   double x, y;
} Vector2;
                     
/* Linear interpolation function. alpha is the position (alpha = 0 takes
   only starting component, alpha = 1 takes only ending component) */
inline double linearInterpolate(double fstart, double fend, double alpha)
{
    return fstart + alpha * (fend - fstart);
}

/* Smoothing function. Takes an input co-ordinate and smoothes it, in the
 * sense that proximity to 0 or 1 pushes the value closer to that extreme
 * limit. The function used is:
     x_smooth = 6 x^5 - 15 x^4 + 10 x^3                  */
inline double smoothInterpolate(double fstart, double fend, double alpha)
{
    double alpha_smoothed = alpha * alpha * alpha * (10 - alpha * (15 - 6 * alpha) );
    return linearInterpolate(fstart, fend, alpha_smoothed);
}

/* Dot product function. Takes the dot product of two 2D vectors */
inline double dotProduct2D(Vector2 v1, Vector2 v2)
{
    return (v1.x * v2.x + v1.y * v2.y);
}

/* Random vector selectio nfunction. Takes as input two integers, and uses
   hashing to convert these into a 'random' selection from 0-255 */
inline int selectVector(unsigned int x, unsigned int y, int *P)
{
    
    //printf("%d\n", P[(x % 256)^P[(y % 256)]]);
    //fflush(stdout);
    return P[(x % 256)^P[(y % 256)]];
}

/* Noise calculation function. Takes an input point in 2D space and
 calculates the dot product with the four vectors surrounding this point */
double noise2D(double x, double y, int *P)
{

    // Convert the input x and y into the corresponding integers
    int int_x = (int)floor(x);
    int int_y = (int)floor(y);
    // Use these to find the "box co-ordinates (co-ords relative to bottom corner of box)
    double box_x = x - int_x;
    double box_y = y - int_y;
    
    // Define vectors to the corners
    Vector2 v_dl = {-box_x, -box_y};
    Vector2 v_dr = {1-box_x, -box_y};
    Vector2 v_ul = {-box_x, 1-box_y};
    Vector2 v_ur = {1-box_x, 1-box_y};
    
    // Use hash function to assign gradient vectors to the corners
    unsigned int k_dl = selectVector(int_x, int_y, P);
    unsigned int k_dr = selectVector(int_x+1, int_y, P);
    unsigned int k_ul = selectVector(int_x, int_y+1, P);
    unsigned int k_ur = selectVector(int_x+1, int_y+1, P);
    // INTEGERS ARE NEGATIVE
    //printf("%d",int_y); 
    //if (int_x <= 0) {
	//char *err = "xxx";
	//printf("%s", err);
        //fflush(stdout);
    //}
    // y definitly prints, and with a retun 0 we can complete the program 
    // so the issue is that I couldn't print k_dl
    // int_x, int_y, int_P exist, so what's happeningg in select vector???
    //fflush(stdout);
    Vector2 g_dl = {vecs_x[k_dl], vecs_y[k_dl]}; 
    Vector2 g_dr = {vecs_x[k_dr], vecs_y[k_dr]};
    Vector2 g_ul = {vecs_x[k_ul], vecs_y[k_ul]};
    Vector2 g_ur = {vecs_x[k_ur], vecs_y[k_ur]};
    
    // Calculate the values of the dot products
    double dl_val = dotProduct2D(v_dl, g_dl);
    double dr_val = dotProduct2D(v_dr, g_dr);
    double ul_val = dotProduct2D(v_ul, g_ul);
    double ur_val = dotProduct2D(v_ur, g_ur);

    // Perform linear interpolation of these dot products and return the result
    double bottom_noise = smoothInterpolate(dl_val, dr_val, box_x);
    double top_noise = smoothInterpolate(ul_val, ur_val, box_x);
    return smoothInterpolate(bottom_noise, top_noise, box_y);
    
}


/* Main runner function. Takes as input the points where Perlin noise is to
 * be evaluated, and calculates the value of octave noise at each point,
 * which is the combination of multiple 'octaves' of Perlin noise. The
 * result is scaled back to lie within [0,1] according to the theoretical
 * maximum and minimum values of 2D ocsity, N_patterns, mesh);
  File "/h1/bmiller/Documents/Perlin_test/lib/genetave noise. */
void OctaveNoise2D(double *points, double *noisefield, unsigned int n, unsigned int N_freq, double fade_factor, int *Ps, double *offsets)
{
    
    // Initialise loop variables
    int i;
    unsigned int j;
    
    // Initialise multipliers
    int freq_mult;
    double scale_mult;
    
    // Initialise individual permutation tables
    int *P;
    
    // Calculate the scaling factor for normalisation
    // This is calculated using a geometric progression, with starting value sqrt(2)/2
    double normalisation_factor = 0.70710678 * (1 - (double)pow((double)fade_factor, (int)N_freq) ) / (1 - fade_factor);
    double normalisation_denominator = 0.5 / normalisation_factor;
    
    // Perform a loop over all input points
    for (i=0; i<n; i++) {
        
        // Initialise value of noisefield here as zero
        noisefield[i] = 0;
        
        // Loop over frequencies
        for (j=0; j<N_freq; j++) {
         
            // Set up multiplier for this frequency
            freq_mult = (int) pow((double)2,(int)j);
            scale_mult = (double)pow((double)fade_factor,(int)j);
            
            // Grab out this octave's permutation table (pointer pointing to shifted position along the full array Ps)
            P = &Ps[256*j];
                        
            // Add value of noisefield at this location
            noisefield[i] = noisefield[i] + scale_mult*noise2D(points[2*i]*freq_mult - offsets[2*j], points[2*i+1]*freq_mult - offsets[2*j+1], P);
            
        }
        
        // Convert this value of the noisefield to its [0,1] normalised equivalent
        noisefield[i] = ( noisefield[i] + normalisation_factor ) * normalisation_denominator;
        
    }
}


// args is a pointer to a python tuple object with each item being an argument 
static PyObject * Octave2D(PyObject * self, PyObject * args)
{

	// instantiate intermediate pointers 
	PyObject * pointsin;
        PyObject * pointsin2;	
	double * points;
	int n_points;
	int i;
	int n_freqs;
	unsigned int N_freq;
	double fade_factor;
	PyObject * Psin;
	int * Ps; 
	PyObject * offsetsin;
	double * offsets;


  	// parse arguments
  	if (!PyArg_ParseTuple(args, "OIdOO", &pointsin2, &N_freq, &fade_factor, &Psin, &offsetsin)) {
    		return NULL;
  	} 
	pointsin = PySequence_Fast(pointsin2, "arguments must be iterable");
	if(!pointsin){
		return 0;
	}
	offsetsin = PySequence_Fast(offsetsin, "arguments must be iterable");
	if(!offsetsin) {
		return 0;
	}
	Psin = PySequence_Fast(Psin, "arguments must be iterable");
	if(!Psin) {
		return 0;
	}
	Py_DECREF(pointsin2);
	
	// prepare points data to array
	n_points = PySequence_Fast_GET_SIZE(pointsin)/2;
	points = malloc(n_points*2*sizeof(double));
        if(!points){
		Py_DECREF(pointsin);
		return PyErr_NoMemory( );
	}	
	for (i=0; i<n_points*2; i++) {
		PyObject *fitem;
		PyObject *item = PySequence_Fast_GET_ITEM(pointsin, i);
		if(!item) {
			Py_DECREF(pointsin);
			free(points);
			return 0;
		}
		fitem = PyNumber_Float(item);
		if(!fitem) {
			Py_DECREF(pointsin);
			free(points);
			PyErr_SetString(PyExc_TypeError, "all items must be numbers");
			return 0;
		}
		points[i] = PyFloat_AS_DOUBLE(fitem);
		Py_DECREF(fitem);
	}
	Py_DECREF(pointsin);

	// prepare offsets data to array 
	n_freqs = PySequence_Fast_GET_SIZE(offsetsin)/2;
	offsets = malloc(n_freqs*2*sizeof(double));
	if(!offsets) {
		Py_DECREF(offsetsin);
		return PyErr_NoMemory( );
	}
	for (i=0; i<n_freqs*2; i++) {
		PyObject *fitem;
		PyObject *item = PySequence_Fast_GET_ITEM(offsetsin, i);
		if(!item) {
			Py_DECREF(offsetsin);
			free(offsets);
			return 0;
		}
		fitem = PyNumber_Float(item);
		if(!fitem) {
			Py_DECREF(offsetsin);
			free(offsets);
			PyErr_SetString(PyExc_TypeError, "all items must be numbers");
			return 0;
		}
		offsets[i] = PyFloat_AS_DOUBLE(fitem);
		Py_DECREF(fitem);
	}
	Py_DECREF(offsetsin);

	// prepare points Ps to array 
	Ps = malloc(n_freqs*256*sizeof(int));
	if(!Ps) {
		Py_DECREF(Psin);
		return PyErr_NoMemory( );
	}
	for (i=0; i<n_freqs*256; i++) {
		PyObject *fitem;
		PyObject *item = PySequence_Fast_GET_ITEM(Psin, i);
		if(!item) {
			Py_DECREF(Psin);
			free(Ps);
			return 0;
		}
		fitem = PyNumber_Float(item);
		if(!fitem) {
			Py_DECREF(Psin);
			free(Ps);
			PyErr_SetString(PyExc_TypeError, "all items must be numbers");
			return 0;
		}
		Ps[i] = PyLong_AsLong(fitem);
		Py_DECREF(fitem);
	}
	Py_DECREF(Psin);

  	// run the actual function 
  	double noisefield[n_points];
	//printf("\n\nmain func runs\n");
	/*
	printf("\n\npoints:\n");
	for (i=0;i<20;i++){
		printf("%f  ", points[i]);
	}	
        printf("\n\nPs:\n");
	for (i=0;i<20;i++){
		printf("%d  ", Ps[i]);
	}	
	printf("\n\noffsets:\n");
	for (i=0;i<10;i++){
		printf("%f  ", offsets[i]);
	}	
	fflush(stdout);
	*/
  	OctaveNoise2D(points, noisefield, n_points, N_freq, fade_factor, Ps, offsets);
	/*
	printf("\n\nnoisefield:\n");
	for (i=0;i<20;i++){
		printf("%f  ", noisefield[i]);
	}	
	fflush(stdout);
	*/
	free(points);
	free(Ps);
	free(offsets);

	// return
	PyObject * lst = PyList_New(n_points);
	if (!lst) {
		return NULL;
	}
	//printf("\n\nnoisefield NUM:\n");
	for (i = 0; i< n_points; i++) {
		PyObject *num = PyFloat_FromDouble(noisefield[i]);
		//if (i<20) {
		//	printf("%f  ", num);
		//	fflush(stdout);
		//}
		if (!num) {
			Py_DECREF(lst);
			return NULL;
		}
		PyList_SET_ITEM(lst, i, num);
	}
  	return lst; 
	//Py_BuildValue("f", lst);
}

static PyMethodDef Octave2DMethods[] = {
	// printf("methods launch");
	{"Octave2D", Octave2D, METH_VARARGS, "Generate Perlin Noise."},
	{NULL, NULL, 0, NULL} /* sentinel */
};

static struct PyModuleDef moduledef = {
	PyModuleDef_HEAD_INIT,
	"Octave3D",
	NULL,
	-1,
	Octave2DMethods
};

//void
PyMODINIT_FUNC PyInit_Octave2D(void)
{
	//printf("init launch");
	//(void) Py_InitModule("Octave2D", Octave2DMethods);
	return PyModule_Create(&moduledef);
}
