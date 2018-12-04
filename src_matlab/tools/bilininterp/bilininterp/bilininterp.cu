#include "mex.h"
#include "cuComplex.h"
#include "bilininterp_kernel.cu"

#define DEBUG_MODE 0 // 1 or 0
#include "bilininterp_defines.h"

// Custom Settings
#define GPU_CUDA_MAX_THREADS_PER_BLOCK 512
typedef unsigned int uint;


//////////////////////////////////////////////////////////////////////////////////////
// helper functions
//////////////////////////////////////////////////////////////////////////////////////
inline uint divUp(uint a, uint b) { return (a + b - 1) / b; }

//////////////////////////////////////////////////////////////////////////////////////
// Forward declarations
//////////////////////////////////////////////////////////////////////////////////////
bool initGPU();
template <typename Type>
void interpolate2d(Type *input, uint input_width, uint input_height, float2 *output_positions, uint output_width, uint output_height, Type *output);
void initCustomMatlabOutputPositions(float2 *output_positions, const mxArray *xi, const mxArray *yi);

//////////////////////////////////////////////////////////////////////////////////////
// MEX Function
//////////////////////////////////////////////////////////////////////////////////////
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
	// init GPU
	if (!initGPU()) {
		mexErrMsgTxt("No CUDA enabled GPU found");
	}

	const mxArray *argZ, *argXI, *argYI;

	// Input checks
	if (!mxIsDouble(prhs[0])) {
		mexErrMsgTxt("Type of input parameter Z must be double, numeric or complex");
	}
	argZ = prhs[0];
	if (!mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2])) {
		mexErrMsgTxt("Type of input parameter XI and YI must be double");
	}
	argXI = prhs[1];
	argYI = prhs[2];
			
	// Output parameter check
	if (nlhs != 1) {
		mexErrMsgTxt("Must have one output argument");
	}

	uint input_width = mxGetN(argZ);
	uint input_height = mxGetM(argZ);
	uint num_input = input_width * input_height;

	uint output_width = mxGetN(argXI);
	uint output_height = mxGetM(argYI);
	uint num_output = output_width * output_height;

	// positions of input points
	float2 *output_positions;
	CUDA_SAFE_CALL( cudaMallocHost((void**)&output_positions, num_output * sizeof(float2)) );
	initCustomMatlabOutputPositions(output_positions, argXI, argYI);

	if (mxIsComplex(argZ))
	{
		//init interpolation output for matlab
		plhs[0] = mxCreateDoubleMatrix(output_height, output_width, mxCOMPLEX);
		double *output_real = mxGetPr(plhs[0]);
		double *output_imag = mxGetPi(plhs[0]);
		if (num_output > 0)
		{
			double *input_real = mxGetPr(argZ);
			double *input_imag = mxGetPi(argZ);
			cuComplex *tmp_input;
			CUDA_SAFE_CALL( cudaMallocHost((void**)&tmp_input, num_input * sizeof(cuComplex)) );
			//convert double input to float/complex
			for (uint y=input_height; y--; ) {
				for (uint x=input_width; x--; ) {
					tmp_input[y * input_width + x] = make_cuComplex((float)input_real[x * input_height + y],
						(float)input_imag[x * input_height + y]);
				}
			}
			
			//reserve space for gpu interpolation output
			cuComplex *tmp_output;
			CUDA_SAFE_CALL( cudaMallocHost((void**)&tmp_output, num_output * sizeof(cuComplex)) );

			//interpolate
			interpolate2d(tmp_input, input_width, input_height, output_positions, output_width, output_height, tmp_output);

			// Convert to double, and we're done
			for (uint y=output_height; y--; ) {
				for (uint x=output_width; x--; ) {
					output_real[x * output_height + y] = (double)cuCrealf(tmp_output[y * output_width + x]);
					output_imag[x * output_height + y] = (double)cuCimagf(tmp_output[y * output_width + x]);
				}
			}

			// Debug Output
			DEBUG_PRINT_CUCOMPLEX_FIELD("Input(real):", tmp_input, input_width, input_height);
			DEBUG_PRINT_FLOAT2_FIELD("Positions:", output_positions, output_width, output_height);
			DEBUG_PRINT_CUCOMPLEX_FIELD("Result:", tmp_output, output_width, output_height);

			// clean up
			CUDA_SAFE_CALL( cudaFreeHost(tmp_input) );
			CUDA_SAFE_CALL( cudaFreeHost(tmp_output) );
		}
	} else {
		//init interpolation output for matlab
		plhs[0] = mxCreateDoubleMatrix(output_height, output_width, mxREAL);
		double *output = mxGetPr(plhs[0]);

		if (num_output > 0)
		{
			//convert double input to float
			double *input = mxGetPr(argZ);
			float *tmp_input;
			CUDA_SAFE_CALL( cudaMallocHost((void**)&tmp_input, num_input * sizeof(float)) );
			for (uint y=input_height; y--; ) {
				for (uint x=input_width; x--; ) {
					tmp_input[y * input_width + x] = (float)input[x * input_height + y];
				}
			}
			
			//reserve space for gpu interpolation output
			float *tmp_output;
			CUDA_SAFE_CALL( cudaMallocHost((void**)&tmp_output, num_output * sizeof(float)) );

			//interpolate
			interpolate2d(tmp_input, input_width, input_height, output_positions, output_width, output_height, tmp_output);

			// Convert back to double
			for (uint y=output_height; y--; ) {
				for (uint x=output_width; x--; ) {
					output[x * output_height + y] = (double)tmp_output[y * output_width + x];
				}
			}

			// Debug Output
			DEBUG_PRINT_FLOAT_FIELD("Input:", tmp_input, input_width, input_height);
			DEBUG_PRINT_FLOAT2_FIELD("Positions:", output_positions, output_width, output_height);
			DEBUG_PRINT_FLOAT_FIELD("Result:", tmp_output, output_width, output_height);
			
			// clean up
			CUDA_SAFE_CALL( cudaFreeHost(tmp_input) );
			CUDA_SAFE_CALL( cudaFreeHost(tmp_output) );
		}
	}
	// clean up
	cudaFreeHost(output_positions);
}


//////////////////////////////////////////////////////////////////////////////////////
// Select a CUDA-ENABLED GPU
// (if there's more than one gpu, the one with the most multiprocessors will be choosen)
//////////////////////////////////////////////////////////////////////////////////////
bool initGPU()
{
	bool result = false;
	int dCnt = 0;
	int selectedCudaDeviceId = 0;
	CUDA_SAFE_CALL( cudaGetDeviceCount(&dCnt) );
	DEBUG("number of cuda gpu devices: %d\n", dCnt);
	if (dCnt > 0) {
		if (dCnt > 1) {
			int multiprocessor_cnt = 0;
			cudaDeviceProp prop;
			for (int deviceId=0; deviceId<dCnt; ++deviceId) {
				if (cudaSuccess == cudaGetDeviceProperties(&prop, deviceId)) {
					if (prop.multiProcessorCount > multiprocessor_cnt) {
						multiprocessor_cnt = prop.multiProcessorCount;
						selectedCudaDeviceId = deviceId;
					}
				}
			}
		} else {
			selectedCudaDeviceId = 0;
		}
		DEBUG("selected device with most multiprocessors: %d\n", selectedCudaDeviceId);
		cudaSetDevice(selectedCudaDeviceId);
		result = true;
	}
	return result;
}


//////////////////////////////////////////////////////////////////////////////////////
void initCustomMatlabOutputPositions(float2 *output_positions, const mxArray *xi, const mxArray *yi)
{
	uint output_width = mxGetN(xi);
	uint output_height = mxGetM(yi);
	double *dxi = mxGetPr(xi);
	double *dyi = mxGetPr(yi);

	uint idx, idx2;
	for (uint x=output_width; x--; ) {
		for (uint y=output_height; y--; ) {
			idx = y * output_width + x;
			idx2 = x * output_height + y;
			output_positions[idx].x = (float)dxi[idx2] - 0.5f;
			output_positions[idx].y = (float)dyi[idx2] - 0.5f;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////
// TEMPLATED INTERPOLATION
//////////////////////////////////////////////////////////////////////////////////////
template <typename Type>
void interpolate2d(Type *input, uint input_width, uint input_height, float2 *output_positions, uint output_width, uint output_height, Type *output)
{
	uint num_output = output_width * output_height;

	cudaArray *dev_input;
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<Type>();
	
	CUDA_SAFE_CALL( cudaMallocArray(&dev_input, &channelDesc, input_width, input_height) );
	CUDA_SAFE_CALL( cudaMemcpyToArray(dev_input, 0, 0, input, input_width * input_height * sizeof(Type), cudaMemcpyHostToDevice) );

	texture<Type, 2, cudaReadModeElementType> &bilinTexture = getTexture<Type>();

	bilinTexture.filterMode = cudaFilterModeLinear;
	bilinTexture.normalized = false;
	CUDA_SAFE_CALL( cudaBindTextureToArray(bilinTexture, dev_input) );

	float2 *dev_output_positions;
	CUDA_SAFE_CALL( cudaMalloc((void**)&dev_output_positions, num_output * sizeof(float2)) );
	CUDA_SAFE_CALL( cudaMemcpy(dev_output_positions, output_positions, num_output * sizeof(float2), cudaMemcpyHostToDevice) );

	Type *dev_output;
	CUDA_SAFE_CALL( cudaMalloc((void**)&dev_output, num_output * sizeof(Type)) );

	// dim3 dimBlock(16, 16);
	//dim3 dimGrid(divUp(output_width, dimBlock.x), divUp(output_height, dimBlock.y));
	dim3 dimBlock(GPU_CUDA_MAX_THREADS_PER_BLOCK, 1);
	dim3 dimGrid(divUp(num_output, dimBlock.x), 1);

	// Interpolation Kernel
	bilinKernel<Type><<<dimGrid, dimBlock>>>(dev_output_positions, dev_output, num_output);

	CUDA_SAFE_CALL( cudaMemcpy(output, dev_output, num_output * sizeof(Type), cudaMemcpyDeviceToHost) );
	CUDA_SAFE_CALL( cudaFree(dev_output) );

	CUDA_SAFE_CALL( cudaUnbindTexture(bilinTexture) );
	CUDA_SAFE_CALL( cudaFreeArray(dev_input) );

	CUDA_SAFE_CALL( cudaFree(dev_output_positions) );
}
