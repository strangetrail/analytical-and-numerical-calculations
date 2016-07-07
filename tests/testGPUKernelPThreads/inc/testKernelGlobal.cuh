//# vi:syntax=cuda

#ifndef _TESTKERNELGLOBAL_CUH_
#define _TESTKERNELGLOBAL_CUH_

// TODO : move all consts to const memory.
// (OPTIONAL. CRITICAL ONLY WITH FDTD IMPLEMENTED INSIDE KERNEL)

// TODO (DONE) : FIX ANOTHER ERROR!!!                                   \
//                Shared memory is different for 2 different streams!!!  
extern __shared__ char memPack [];

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

  do
  {
    for ( k = 0; k < args.maxChunks; k++ )
    {
      for ( l = 0; l < dimzBlock; l++ )
      {
        for ( i = 0; i < dimxBlock; i++ )
          for ( j = 0; j < dimyBlock; j++ )
            do
            {
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

                  iGRFdevice = deviceGlobalRefreshFlags[idx_io];
                  i *= iGRFdevice;// TODO : ` i = i & ( (~f ^ f) | f ); '
                  j *= iGRFdevice;
                }
              }
            }
            while ( !iGRFdevice );

        // `iz' or `l' - is a XxY blocks slice index in global memory:
        hostWait4RefreshGlobalSlice[l] = 0;
      }
      // TODO : This is much more optimal than overflow technique,   \
                 because it has lesser number of conditional checks:  
      for ( l = 0; l < dimzBlock; l++ )
        l *= hostWait4RefreshGlobalSlice[l];
      *hostWait4RefreshingChunk_WhileLoadingSlices = 0;
      // TODO : Do I need to block acces to volatile host-device memory \
      //         section?? ANSWER : AT LEAST YOU NEED VOLATILE          \
      //         MODIFICATOR:                                            
      while ( *hostWaitWhileLoadingGlobalChunk ) {}
      *hostWaitWhileLoadingGlobalChunk = 1;
      *deviceWaitWhileLoadingGlobalChunk = 0;
    }
  }
  while ( *( args.bContinue ) );
}

__global__ void testKernel ( TestKernelArguments_t args )
{
  /*
  extern __shared__ float tile [];
  */

  // TODO : Remove extra `const':
  const int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
  // TODO : Optimize!!!
  const int gtidy = blockIdx.y * blockDim.y + threadIdx.y;
  // TODO : Optimize!!!
  const int ltidx = threadIdx.x;
  const int ltidy = threadIdx.y;
  const int blkx = blockIdx.x;
  const int blky = blockIdx.y;
  unsigned char *deviceWaitWhileLoadingGlobalChunk =      \
                  args.deviceWaitWhileLoadingGlobalChunk;  
  const int /*&dimx      = args.dimx,*/
            /* dimx in whole XxY blocks global memory slices !!! */
            /* e.g. slice No 1, slice No 2, ... , slice No N */
            /*&dimy      = args.dimy,*/
            /* dimy in whole XxY blocks global memory slices !!! */
            &dimz      = args.dimz,
            /* `dimz' is thread memory length in z direction !!! */
            &dimxBlock = args.dimxBlock,
            &dimyBlock = args.dimyBlock,
            &dimzBlock = args.dimSlice,
            &dimThreadsX = args.dimThreadsX,
            &dimThreadsY = args.dimThreadsY;
  // TODO : Does it correct to get by reference from struct passed as a \
  //         parameter to cuda kernel?                                   
  unsigned char *deviceGlobalRefreshFlags = args.deviceGlobalRefreshFlags;

  int iz, iiz, idx_io;

  // TODO : ALIGN CODE LINES!!!
  float  fResult,
        *ioBuffer = args.buffer,
        /* ioBuffer[dimSlice][dimxBlock][dimyBlock] */
        /*         [dimThreadsX][dimThreadsY]       */
        *ioTile   = (float *)memPack;
        /* ioTile[dimThreadsX][dimThreadsY] */

  clock_t    start      = clock(),
             now,
             cycles,
          * &global_now = args.global_now;

  // TODO : Put here global loop over all global memory CHUNCKS \
  //         (index is irrelevant - host responsible for        \
  //         proper memory loading and unloading):               
  do
  {

    // Per-block loop alongside z direction:
#pragma unroll 3
    for ( iz = 0; iz < dimzBlock; iz++ )
    {

      idx_io = iz * dimxBlock * dimyBlock * dimThreadsX * dimThreadsY \
               + blkx * dimyBlock * dimThreadsX * dimThreadsY         \
               + blky * dimThreadsX * dimThreadsY                     \
               + ltidx * dimThreadsY                                  \
               + ltidy;                                                

      memcpy                                           \
      (                                                \
        (float *)&ioTile[ltidx * dimThreadsY + ltidy], \
        (float *)&ioBuffer[idx_io],                    \
        dimz * sizeof(float)                           \
      );                                                

      __syncthreads();

#pragma unroll 3
      for ( iiz = 0; iiz < dimz - 1/* Optimize !! */; iiz++ )
      {
        /* TODO : I AM HERE (7.07.16) : ERROR: Out-fo-range exception in */
        /*                                      cuda kernel!             */
        fResult = ioTile[ltidy + dimThreadsY*ltidx + dimThreadsY*dimThreadsX*iiz] + ioTile[ltidy + dimThreadsY*ltidx + dimThreadsY*dimThreadsX*(iiz + 1)/* Optimize */];
        for (;;)
        {
          now = clock();
          cycles = now > start ? now - start : now + (0xffffffff - start);
          if ( cycles >= 10000 )
          {
            break;
          }
        }
        *global_now = now;
        fResult *= ioTile[ltidx*dimThreadsY + ltidy + (iiz)*dimThreadsX*dimThreadsY] - ioTile[ltidx*dimThreadsY + ltidy + (iiz + 1)*dimThreadsX*dimThreadsY/* Optimize */];
        __syncthreads();
        ioTile[ltidx*dimThreadsY + ltidy + (iiz + 1)*dimThreadsX*dimThreadsY] = fResult;
        __syncthreads();
      }

      memcpy ( (float *)&ioBuffer[idx_io],                    \
               (float *)&ioTile[ltidx * dimThreadsY + ltidy], \
               dimzBlock * sizeof(float)                      \
             );                                                

      // TODO (DONE) : FIX AN ERROR!!! Flags have to include both \
      //                                xy-block and z indexes.    
      // TODO (DONE) : Find out why there are only thread indices but no \
      //                block.                                            
      // TODO : We need control stream for synchronization with host:
      deviceGlobalRefreshFlags[idx_io] = 1;
      // 1 - ready, 0 - wait.

    }

    __syncthreads();

    while ( *deviceWaitWhileLoadingGlobalChunk ) {}
    // ???????????? TODO : how to reset???

    __syncthreads();//???? TODO : remove??? __syncthreads()

    *deviceWaitWhileLoadingGlobalChunk = __any(1);

    // TODO : Verify if this reset works properly.
    for ( iz = 0 ; iz < dimz ; iz++ )
      deviceGlobalRefreshFlags[gtidx * dimThreadsY + gtidy] = 0;
      // 1 - ready, 0 - wait

    //Stop and wait here in GLOBAL loop for new CHUNK \
    // (synchronize with host)                         
  }
  while ( *( args.bContinue ) );
  // TODO : Verify repeat-until behavior in c++.

}

#endif
