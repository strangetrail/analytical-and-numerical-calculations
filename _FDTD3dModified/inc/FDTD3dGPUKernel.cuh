//# vi:syntax=cuda
/*
 * Copyright 1993-2015 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

#ifndef _FDTD3DGPUKERNEL_CUH_
#ifndef _FDTD3DGPUKERNEL_CUH_

#include "FDTD3dGPU.h"

// Note: If you change the RADIUS or VOLUME, you should also change the unrolling below
#define VOLUME VOLUME_SHARED
#define RADIUS RADIUS_SHARED

__constant__ UpdateCoefficients_t stencil;
__constant__ int T_src;
__constant__ float omega_src;
__constant__ float n_inc;
__constant__ float k_src;
__constant__ float musqrt_src;
__constant__ float epsilonsqrt_src;
__constant__ int blockDimXY;

// TODO : Move all consts to const memory.              \
//         CRITICAL WITH FDTD IMPLEMENTED INSIDE KERNEL~~
// TODO (DONE) : CRITICAL : FIX ANOTHER ERROR!!!            \
//                           Shared memory is different for \
//                           2 different streams!!!~~~~~~~~~~
// TODO (DONE) : Threads, pthreads, streams, or processes \
//                hangs.~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
extern __shared__ char memPack []; /* Test: 16 * 16 * 24 bytes. */
extern __shared__ float tile [];

// TODO : Ensure that control stream do not wait for threads while they \
//         are reach chunk end and only waits for slices!!! Send        \
//         signals for slice reloading after each slice finished. Also  \
//         ensure that threads will wait after finishing with chunk     \
//         until new chunk is copletelly loaded by per-slice pthreads.   
__global__ void testKernelControlStream                   \
                (                                         \
                  const TestControlKernelArguments_t args \
                )                                          
{
  const int Xthreads = args.dimThreadsX,
            Ythreads = args.dimThreadsY,
            dimzBlock = args.dimSlice,
            dimyBlock = args.dimyBlock,
            dimxBlock = args.dimxBlock;

  // TODO : You need to use multiplications while indexxing elements.    \
  //         Otherwise you may want to use devceQuery in Makefile before \
  //         compiling code to link it with proper static header with    \
  //         particular GPU architecture description. There is no other  \
  //         way to optimize dynamic array indexing.                      
  // TODO : Comment each flag that people understand what flags are set \
  //         for.                                                        
  unsigned char *deviceWaitWhileLoadingGlobalChunk =      \
                  args.deviceWaitWhileLoadingGlobalChunk,  
                *deviceGlobalRefreshFlags = args.deviceGlobalRefreshFlags,
                *hostWaitWhileLoadingGlobalChunk =      \
                  args.hostWaitWhileLoadingGlobalChunk,  
                *hostWait4RefreshGlobalSlice =      \
                  args.hostWait4RefreshGlobalSlice,  
                *hostWait4RefreshingChunk_WhileLoadingSlices =      \
                  args.hostWait4RefreshingChunk_WhileLoadingSlices;  

  // Peviously there are was a rule that initial values of `i' and `j' \
  //  should be `MAXSHORTINT - 1 - (dimThreadsX | dimThreadsY)';       \
  //  but now it's OBSOLETE.                                            
  int i, j, k, l, m, n, iGRFdevice, idx_io;

  /* STEP4 : Waiting while main process sets continuation flag. */
  // TODO (DONE) : FIX CRITICAL ERROR : Loop does not count \
  //                                     timesteps:          
  // TODO I AM HERE (09.03.16) : FIX AN ERROR : Implement an empty loop     \
  //                                             that waits for `bContinue' \
  //                                             flag:                       
  // TODO : Use loop with timesteps instead of bContinue flag:
  while ( !*(args.bContinue) ) {}
  for ( int timestep = 0; timestep < args.timesteps; timestep++ )
  {
    for ( k = 0; k < args.maxChunks; k++ )
    {
      for ( l = 0; l < dimzBlock; l++ )
      {
        for ( i = 0; i < dimxBlock; i+=iGRFdevice )
          for ( j = 0; j < dimyBlock; j+=iGRFdevice )
            do
            {
              iGRFdevice = 1;
              for ( m = 0; m < Xthreads; m++ )
              {
                for ( n = 0; n < Ythreads; n++ )
                {
                  idx_io = l * dimxBlock * dimyBlock             \
                           * Xthreads * Ythreads                 \
                           + i * dimyBlock * Xthreads * Ythreads \
                           + j * Xthreads * Ythreads             \
                           + m * Ythreads                        \
                           + n;                                   

                  /* STEP9 : Wait until all threads in all blocks sets flags */
                  /*          indicating that next slice can be reloaded     */
                  /*          from host memory to device memory.             */
                  iGRFdevice *= deviceGlobalRefreshFlags[idx_io];
                  // TODO : ` i = i & ( (~f ^ f) | f ); '
                }
              }
            }
            while ( !iGRFdevice );

        /* STEP14 : Control stream sets flag telling main process that it   */
        /*           can copy slice from the device memory to host and load */
        /*           new data into device memory.                           */
        // `iz' or `l' - is a XxY blocks slice index in global memory:
        hostWait4RefreshGlobalSlice[l] = 0;
      }

      /* STEP18 : Waiting until PThreads reload all slices from device to */
      /*           host memory and load new slices to device.             */
      // TODO : This is much more optimal than overflow technique,   \
                 because it has lesser number of conditional checks:  
      // TODO I AM HERE (07.30.16) : Test the loop below:
      l = 0;
      while ( l < dimzBlock )
        l += hostWait4RefreshGlobalSlice[l];

      /* STEP19 : Telling host process that global chunk can be reloaded, */
      /*           evaluated, updated, etc.                               */
      *hostWait4RefreshingChunk_WhileLoadingSlices = 0;

      // TODO : Do I need to block acces to volatile host-device memory \
      //         section?? ANSWER : AT LEAST YOU NEED VOLATILE          \
      //         MODIFICATOR:                                            
      /* STEP21 : Waiting until global chunk is reloaded. */
      while ( *hostWaitWhileLoadingGlobalChunk ) {}

      /* STEP22 : Resetting flag to its "wait" state. */
      *hostWaitWhileLoadingGlobalChunk = 1;

      /* STEP23 : Telling kernel threads to continue calculations. */
      *deviceWaitWhileLoadingGlobalChunk = 0;

      /* STEP28 : Wait for first kernel thread reached the line telling   */
      /*           others that all flags, indicating that specific global */
      /*           memory part needs to be updated in shared memory, have */
      /*           been reset for all particular blocks and threads.      */
      while ( !(*deviceWaitWhileLoadingGlobalChunk) ) {}
    }
  }
}

__global__ void FiniteDifferencesKernelWithSource          \
                (                                          \
                  FDTDKernelWithSourceArguments_t args     \
                  /* TODO : Pass sinsrc with input data */ \
                )                                           
{
    bool validr = true;
    bool validw = true;
    // TODO : Remove extra `const':
    // TODO : Optimize!!!
    const int /*gtidx = blockIdx.x * blockDim.x + threadIdx.x,*/
              /*gtidy = blockIdx.y * blockDim.y + threadIdx.y,*/
              ltidx = threadIdx.x,
              ltidy = threadIdx.y,
              blkx  = blockIdx.x,
              blky  = blockIdx.y;

    // TODO : Does it correct to get values by reference from struct \
    //         that passed as a parameter to cuda kernel?             
    unsigned char *deviceWaitWhileLoadingGlobalChunk =      \
                    args.deviceWaitWhileLoadingGlobalChunk,  
                  *deviceGlobalRefreshFlags = args.deviceGlobalRefreshFlags;
    // TODO (IMPORTANT and INTERESTING) : Figure out why reference types \
    //                                     behaves as register pointers  \
    //                                     in CUDA:                       
    const int /*&dimx          = args.dimx,*/
              /* `dimx' is whole XxY blocks global memory slices !!! */
              /*  e.g. slice No 1, slice No 2, ... , slice No N.     */
              /*&dimy          = args.dimy,*/
              /* dimy in whole XxY blocks global memory slices !!! */
              /*&*/dimz        = args.dimz,
              /* `dimz' is a thread memory length in z direction !!! */
              /*&*/dimxBlock   = args.dimxBlock,
              /*&*/dimyBlock   = args.dimyBlock,
              /* `dimxBlock' and `dimyBlock' is a length of grid in blocks */
              /*  alongside x and y direction respectively.                */
              /*&*/dimzBlock   = args.dimSlice,
              /* `dimzBlock' is a number of slices per chunk. */
              /*&*/dimThreadsX = args.dimThreadsX,
              /*&*/dimThreadsY = args.dimThreadsY;

    int k, iz, iiz, idx_io, idx_sync, idx_shared;

    // TODO : ALIGN CODE LINES!!!
    float  fResult,
          *ioBuffer = args.buffer,
          /* ioBuffer[dimSlice][(*dimxBlock)][(*dimyBlock)] */
          /*         [dimThreadsX][dimThreadsY]             */
          *ioTile   = (float *)memPack;
          /* ioTile[dimThreadsX][dimThreadsY] */

    TFSF_t fdtdTFSFsrcE,
           fdtdTFSFsrcH;

  /* STEP6 : Waiting when main process sets continuation flag. */
  // TODO (DONE) : Put here global loop over all global memory CHUNCKS \
  //                (index is irrelevant - host responsible for        \
  //                proper memory loading and unloading):               
  while ( !*(args.bContinue) ) {}
  for ( int timestep = 0; timestep < args.timesteps; timestep++ )
  {
    for ( k = 0; k < args.maxChunks; k++ )
    {
      // Per-block loop alongside z direction:
//#pragma unroll 3
      for ( iz = 0; iz < dimzBlock; iz++ )
      {
        // ioBuffer size is equal to                               \
        //  dimSlice*gridDim*gridDim*blockDim*blockDim*threadSize.  
        idx_io = dimxBlock*dimyBlock*dimThreadsX*dimThreadsY*dimz*iz \
                 + dimyBlock*dimThreadsX*dimThreadsY*dimz*blkx       \
                 + dimThreadsX*dimThreadsY*dimz*blky                 \
                 + dimThreadsY*dimz*ltidx                            \
                 + dimz*ltidy;                                        
        idx_sync = dimxBlock*dimyBlock*dimThreadsX*dimThreadsY*iz \
                   + dimyBlock*dimThreadsX*dimThreadsY*blkx       \
                   + dimThreadsX*dimThreadsY*blky                 \
                   + dimThreadsY*ltidx                            \
                   + ltidy;                                        

        memcpy                                                               \
        (                                                                    \
          (FieldComponents_t *)&ioTile[dimThreadsY*dimz*ltidx + dimz*ltidy], \
          (FieldComponents_t *)&ioBuffer[idx_io],                            \
          dimz * sizeof(FieldComponents_t)                                   \
        );                                                                    
        /* TODO : I AM HERE (09.10.16) : memcpy ( IH and ID ) likewise */
        /* TODO : memcpy fdtd tf-sf sin sources ??? */

        /* STEP10 : Wait until all threads comes to this line. */
        __syncthreads();

        // Step through xy-planes:
        // Inside PML pseudoboundaries:
        // TF-SF corrections:
        // TODO : Change ix / iy to specific coordinates associated \
        //         with threads and blocks !!!!                      
        fdtdTFSFsrcE.sinsrc = (FieldComponents_t ***)&/*tile???*/;
        fdtdTFSFsrcH.getTFSF = &getFDTDTFSFsrcE;

        fdtdTFSFsrcH.sinsrc = (FieldComponents_t ***)&/*tile???*/[blockDimXY];
        fdtdTFSFsrcH.getTFSF = &getFDTDTFSFsrcH;

        fdtdRef4SingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, dimx, dimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, RADIUS_SHARED, 0, fdtdRefFieldWindowMediaHWrapper, fdtdRefFieldWindowMediaDWrapper, fdtdTFSFsrcE, 0 );

        fdtdRefSingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, RADIUS_SHARED, 0, 0, fdtdRefFieldWindowVolume, fdtdRefFieldWindowMediaDWrapper, fdtdTFSFsrcE, 0 );

        fdtdRef4SingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, dimx, dimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, RADIUS_SHARED + 1, 0, fdtdRefFieldWindowMediaHWrapper, fdtdRefFieldWindowMediaDWrapper, 0, fdtdTFSFsrcH );

        fdtdRefSingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, RADIUS_SHARED + 1, 0, 0, fdtdRefFieldWindowVolume, fdtdRefFieldWindowMediaDWrapper, 0, fdtdTFSFsrcH );

        // Field around sources:
        for ( int iiz = RADIUS_SHARED ; iiz < RADIUS_SHARED + 2; iiz++ )
        {
          fdtdRefSingleXY ( 0, dimXYPostPML, 0, outerDimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, ix, 0, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
          fdtdRefSingleXY ( dimXYPrePML, outerDimx, 0, outerDimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, ix, dimx, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
          fdtdRefSingleXY ( dimXYPostPML, dimXYPrePML, 0, dimXYPostPML, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, iy, 0, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
          fdtdRefSingleXY ( dimXYPostPMLu, dimXYPrePML, dimXYPrePML, outerDimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, iy, dimy, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
        }

        for (int iiz = RADIUS_SHARED + 2 ; iiz < dimZPrePML ; iiz++)
        {
          fdtdRef4SingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, dimx, dimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, 0, fdtdRefFieldWindowMediaHWrapper, fdtdRefFieldWindowMediaDWrapper, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );

          fdtdRefSingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, 0, 0, fdtdRefFieldWindowVolume, fdtdRefFieldWindowMediaDWrapper, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );

          fdtdRefSingleXY ( 0, dimXYPostPML, 0, outerDimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, ix, 0, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
          fdtdRefSingleXY ( dimXYPrePML, outerDimx, 0, outerDimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, ix, dimx, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
          fdtdRefSingleXY ( dimXYPostPML, dimXYPrePML, 0, dimXYPostPML, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, iy, 0, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
          fdtdRefSingleXY ( dimXYPostPMLu, dimXYPrePML, dimXYPrePML, outerDimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, iy, dimy, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
        }

        for (int iiz = 0 ; iiz < RADIUS_SHARED ; iiz++)
        {
          fdtdRef4SingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, dimx, dimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, iiz, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
          fdtdRef4SingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, dimx, dimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, RADIUS_SHARED + dimz + iiz, RADIUS_SHARED - iiz, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );

          fdtdRefSingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, RADIUS_SHARED - iiz, 0, fdtdRefFieldPMLVolumeH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );//RADIUS_SHARED - iiz - that may cause an error due to passing argument inf subrouting by reference and not value --- Need min and max values and not only min value
          fdtdRefSingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz + dimz + RADIUS_SHARED, iiz, 0, fdtdRefFieldPMLVolumeH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
        }
      }
    }
  }
}

__global__ void FiniteDifferencesKernel(float *output,
                                        const float *input,/* place snsrc with input data */
                                        const int dimx,
                                        const int dimy,
                                        const int dimz)
{
    bool validr = true;
    bool validw = true;
    const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
    const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
    const int ltidx = threadIdx.x;
    const int ltidy = threadIdx.y;
    const int workx = blockDim.x;
    const int worky = blockDim.y;
    __shared__ float tile[k_blockDimMaxY + 2 * RADIUS][k_blockDimX + 2 * RADIUS];

    const int stride_y = dimx + 2 * RADIUS;
    const int stride_z = stride_y * (dimy + 2 * RADIUS);

    int inputIndex  = 0;
    int outputIndex = 0;

    // Advance inputIndex to start of inner volume
    inputIndex += RADIUS * stride_y + RADIUS;

    // Advance inputIndex to target element
    inputIndex += gtidy * stride_y + gtidx;

    float infront[RADIUS];
    float behind[RADIUS];
    float current;

    const int tx = ltidx + RADIUS;
    const int ty = ltidy + RADIUS;

    // Check in bounds
    if ((gtidx >= dimx + RADIUS) || (gtidy >= dimy + RADIUS))
        validr = false;

    if ((gtidx >= dimx) || (gtidy >= dimy))
        validw = false;

    // Preload the "infront" and "behind" data
    for (int i = RADIUS - 2 ; i >= 0 ; i--)
    {
        if (validr)
            behind[i] = input[inputIndex];

        inputIndex += stride_z;
    }

    if (validr)
        current = input[inputIndex];

    outputIndex = inputIndex;
    inputIndex += stride_z;

    for (int i = 0 ; i < RADIUS ; i++)
    {
        if (validr)
            infront[i] = input[inputIndex];

        inputIndex += stride_z;
    }

    // Step through the xy-planes

    // Inside PML pseudoboundaries
    // TODO: change ix / iy to specific coordinates associated with threads and blocks !!!!

    //TF-SF corrections

    fdtdTFSFsrcE = 

    fdtdRef4SingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, dimx, dimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, RADIUS_SHARED, 0, fdtdRefFieldWindowMediaHWrapper, fdtdRefFieldWindowMediaDWrapper, &getFDTDTFSFsrcE, &getFDTDTFSFsrcNull );

    fdtdRefSingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, RADIUS_SHARED, 0, 0, fdtdRefFieldWindowVolume, fdtdRefFieldWindowMediaDWrapper, &getFDTDTFSFsrcE, &getFDTDTFSFsrcNull );

    fdtdRef4SingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, dimx, dimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, RADIUS_SHARED + 1, 0, fdtdRefFieldWindowMediaHWrapper, fdtdRefFieldWindowMediaDWrapper, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcH );

    fdtdRefSingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, RADIUS_SHARED + 1, 0, 0, fdtdRefFieldWindowVolume, fdtdRefFieldWindowMediaDWrapper, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcH );

    // Field around sources

    for ( int iiz = RADIUS_SHARED ; iiz < RADIUS_SHARED + 2; iiz++ )
    {
      fdtdRefSingleXY ( 0, dimXYPostPML, 0, outerDimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, ix, 0, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
      fdtdRefSingleXY ( dimXYPrePML, outerDimx, 0, outerDimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, ix, dimx, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
      fdtdRefSingleXY ( dimXYPostPML, dimXYPrePML, 0, dimXYPostPML, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, iy, 0, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
      fdtdRefSingleXY ( dimXYPostPMLu, dimXYPrePML, dimXYPrePML, outerDimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, iy, dimy, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
    }

    for (int iiz = RADIUS_SHARED + 2 ; iiz < dimZPrePML ; iiz++)
    {
      fdtdRef4SingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, dimx, dimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, 0, fdtdRefFieldWindowMediaHWrapper, fdtdRefFieldWindowMediaDWrapper, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );

      fdtdRefSingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, 0, 0, fdtdRefFieldWindowVolume, fdtdRefFieldWindowMediaDWrapper, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );

      fdtdRefSingleXY ( 0, dimXYPostPML, 0, outerDimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, ix, 0, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
      fdtdRefSingleXY ( dimXYPrePML, outerDimx, 0, outerDimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, ix, dimx, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
      fdtdRefSingleXY ( dimXYPostPML, dimXYPrePML, 0, dimXYPostPML, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, iy, 0, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
      fdtdRefSingleXY ( dimXYPostPMLu, dimXYPrePML, dimXYPrePML, outerDimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, iy, dimy, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
    }

    for (int iiz = 0 ; iiz < RADIUS_SHARED ; iiz++)
    {
      fdtdRef4SingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, dimx, dimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, iiz, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
      fdtdRef4SingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, dimx, dimy, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, RADIUS_SHARED + dimz + iiz, RADIUS_SHARED - iiz, fdtdRefFieldPMLMediaH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );

      fdtdRefSingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz, RADIUS_SHARED - iiz, 0, fdtdRefFieldPMLVolumeH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );//RADIUS_SHARED - iiz - that may cause an error due to passing argument inf subrouting by reference and not value --- Need min and max values and not only min value
      fdtdRefSingleXY ( dimXYHalfNearVolume, dimXYHalfAfterVolume, dimXYHalfNearVolume, dimXYHalfAfterVolume, ix, iy, *output, *input, updateCoeffs, IH, ID, Curl, iiz + dimz + RADIUS_SHARED, iiz, 0, fdtdRefFieldPMLVolumeH, fdtdRefFieldPMLMediaD, &getFDTDTFSFsrcNull, &getFDTDTFSFsrcNull );
    }

#pragma unroll 9

    for (int iiz = 0 ; iiz < dimz ; iiz++)
    {
        // Advance the slice (move the thread-front)
        for (int i = RADIUS - 1 ; i > 0 ; i--)
            behind[i] = behind[i - 1];

        behind[0] = current;
        current = infront[0];
#pragma unroll 4

        for (int i = 0 ; i < RADIUS - 1 ; i++)
            infront[i] = infront[i + 1];

        if (validr)
            infront[RADIUS - 1] = input[inputIndex];

        inputIndex  += stride_z;
        outputIndex += stride_z;
        __syncthreads();

        // Note that for the work items on the boundary of the problem, the
        // supplied index when reading the halo (below) may wrap to the
        // previous/next row or even the previous/next xy-plane. This is
        // acceptable since a) we disable the output write for these work
        // items and b) there is at least one xy-plane before/after the
        // current plane, so the access will be within bounds.

        // Update the data slice in the local tile
        // Halo above & below
        if (ltidy < RADIUS)
        {
            tile[ltidy][tx]                  = input[outputIndex - RADIUS * stride_y];
            tile[ltidy + worky + RADIUS][tx] = input[outputIndex + worky * stride_y];
        }

        // Halo left & right
        if (ltidx < RADIUS)
        {
            tile[ty][ltidx]                  = input[outputIndex - RADIUS];
            tile[ty][ltidx + workx + RADIUS] = input[outputIndex + workx];
        }

        tile[ty][tx] = current;
        __syncthreads();

        // Compute the output value
        float value = stencil[0] * current;
#pragma unroll 4

        for (int i = 1 ; i <= RADIUS ; i++)
        {
            value += stencil[i] * (infront[i-1] + behind[i-1] + tile[ty - i][tx] + tile[ty + i][tx] + tile[ty][tx - i] + tile[ty][tx + i]);
        }

        // Store the output value
        if (validw)
            output[outputIndex] = value;
    }
}
