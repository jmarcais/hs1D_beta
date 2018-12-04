#ifndef _BILININTERP_KERNEL_CU_
#define _BILININTERP_KERNEL_CU_

#include "cuComplex.h"

//#define uIntMul(a,b) a * b
#define uIntMul(a,b) __umul24(a, b)

//////////////////////////////////////////////////////////////////////////////////////
// Currently, CUDA requires all textures to be file-scoped
// (@see CUDA Programming Guide 4.3.4.1)
//////////////////////////////////////////////////////////////////////////////////////
texture<float, 2, cudaReadModeElementType> bilinTexFloat;
texture<cuComplex, 2, cudaReadModeElementType> bilinTexComplex;

template <typename Type>
texture<Type, 2, cudaReadModeElementType> &getTexture();

template <typename Type>
static inline __device__ texture<Type, 2, cudaReadModeElementType> &getDeviceTexture();

template <>
texture<float, 2, cudaReadModeElementType> &getTexture<float>() { return bilinTexFloat; }

template <>
static __inline__ __device__ texture<float, 2, cudaReadModeElementType> &getDeviceTexture<float>() { return bilinTexFloat; }

template <>
texture<cuComplex, 2, cudaReadModeElementType> &getTexture<cuComplex>() { return bilinTexComplex; }

template <>
static __inline__ __device__ texture<cuComplex, 2, cudaReadModeElementType> &getDeviceTexture<cuComplex>() { return bilinTexComplex; }


//////////////////////////////////////////////////////////////////////////////////////
// Cuda kernel for interpolation templated
//////////////////////////////////////////////////////////////////////////////////////
template <typename Type>
__global__ void bilinKernel(float2* positions, Type* output, unsigned int num_output)
{
	const int idx = uIntMul(blockIdx.x, blockDim.x) + threadIdx.x;
	if (idx < num_output) {
		output[idx] = tex2D(getDeviceTexture<Type>(), positions[idx].x, positions[idx].y);
	}
}

#endif //#ifndef _BILININTERP_KERNEL_CU_