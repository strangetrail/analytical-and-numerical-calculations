/*
 *  * example.c: very simple example of port I/O
 *   *
 *    * This code does nothing useful, just a port write, a pause,
 *     * and a port read. Compile with `gcc -O2 -o example example.c',
 *      * and run as root with `./example'.
 *       */

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/io.h>

#define BASEPORT 0x3bc /* usblp0 */
/*#define BASEPORT 0x378*/ /* usblp1 */
/*#define BASEPORT 0x278*/ /* usblp2 */

int main()
{
  unsigned char byte;
  short i;
  int count;
  /* Get access to the ports */
  if (ioperm(BASEPORT, 3, 1)) {perror("ioperm"); exit(1);}
  
  /* Set the data signals (D0-7) of the port to all low (0) */
  count = 0;
  byte = 0x00;
  do
  {
    outb(byte, BASEPORT+2);
    outb(byte, BASEPORT);
    /*i = 0;
    for (i = 0; i < 30; i = i + 1)
    {
      usleep(100000);
    }*/
    byte = ~byte;
    count = count + 1;
  }
  while(count < 20000000);
  
  /* Sleep for a while (100 ms) */
  usleep(100000);
  
  /* Read from the status port (BASE+1) and display the result */
  printf("status: %d\n", inb(BASEPORT + 1));

  /* We don't need the ports anymore */
  if (ioperm(BASEPORT, 3, 0)) {perror("ioperm"); exit(1);}

  exit(0);
}

/* end of example.c */
