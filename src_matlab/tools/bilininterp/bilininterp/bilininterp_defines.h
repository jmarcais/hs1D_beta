#ifndef _BILININTERP_DEFINES_H_
#define _BILININTERP_DEFINES_H_

#if DEBUG_MODE
	#define DEBUG(s, ...) mexPrintf(s, __VA_ARGS__) //Variadic Macro
	#define DEBUG_PRINT_FLOAT_FIELD(title, data, width, height) mexPrintf("%s\n",title); for (uint y=0; y<height; ++y) { for (uint x=0; x<width; ++x) { mexPrintf(" [%5.2f]", data[y*width+x]); } mexPrintf("\n"); }
	#define DEBUG_PRINT_FLOAT2_FIELD(title, data, width, height) mexPrintf("%s\n",title); for (uint y=0; y<height; ++y) { for (uint x=0; x<width; ++x) { mexPrintf(" [%5.2f, %5.2f]", data[y*width+x].x, data[y*width+x].y); } mexPrintf("\n"); }
	#define DEBUG_PRINT_CUCOMPLEX_FIELD(title, data, width, height) DEBUG_PRINT_FLOAT2_FIELD(title, data, width, height)
	#define CUDA_SAFE_CALL(call) do {                              \
		cudaError err = call;                                      \
		if( cudaSuccess != err) {                                  \
			DEBUG("Cuda error in file '%s' in line %i : %s.\n",    \
				__FILE__, __LINE__, cudaGetErrorString( err) );    \
		} } while (0)
#else
	#define DEBUG(s, ...)
	#define DEBUG_PRINT_FLOAT_FIELD(title, data, width, height)
	#define DEBUG_PRINT_FLOAT2_FIELD(title, data, width, height)
	#define DEBUG_PRINT_CUCOMPLEX_FIELD(title, data, width, height) DEBUG_PRINT_FLOAT2_FIELD(title, data, width, height)
	#define CUDA_SAFE_CALL(call) call
#endif

#endif //_BILININTERP_DEFINES_H_